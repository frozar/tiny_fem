# -*-coding: utf-8 -*-
import os
import string
import numpy as np
# import pandas as pd
import pylab as pl
import scipy.interpolate

class sample():
    """ Classe définissant un échantillon par :
    - son nom
    - son matériau
    - son maillage """
    
        
    def __init__(self, nom, materiau, maillage):
        "Initialisation de la classe"
        self.name = nom
        self.materiau = materiau
        self.maillage = maillage
        
    def mesh(self, lx, ly, nx, ny):
		
        " CREATION D'UN MAILLAGE RECTANGULAIRE GMSH "
	
	try : 
		os.mkdir('MESH')
	except :
		pass
		
		
	self.geo_file = "MESH" + os.sep + self.name +'.geo'
	geo_file = open(self.geo_file, 'w')

	geo_file.write('/* Description de la geometrie */\n')
	geo_file.write('/* LARGEUR */\n')	
	geo_file.write('largeur = '+str(lx)+';\n')
	geo_file.write('/* LONGUEUR */\n')	
	geo_file.write('longueur = '+str(ly)+';\n')		
	geo_file.write('/* NOMBRE DE NOEUDS SUIVANT X */\n')	
	geo_file.write('nx = '+str(nx+1)+';\n')	
	geo_file.write('/* NOMBRE DE NOEUDS SUIVANT Y */\n')	
	geo_file.write('ny = '+str(ny+1)+';\n')	
	geo_file.write('/* FIN Description de la geometrie */\n\n')		
		
		
	## CREATION DES POINTS DU RECTANGLE ##
	geo_file.write('Point(1) = {0. , 0. , 0. , 1. };\n')
	geo_file.write('Point(2) = {largeur , 0. , 0. , 1. };\n')
	geo_file.write('Point(3) = {largeur , longueur , 0., 1. };\n')
	geo_file.write('Point(4) = {0. , longueur , 0. , 1. };\n')
		
		
	#~ ## CREATION DES LIGNES DU RECTANGLE ##
	geo_file.write('Line(10) = {1, 2}; \n')
	geo_file.write('Line(20) = {2, 3}; \n')
	geo_file.write('Line(30) = {3, 4}; \n')
	geo_file.write('Line(40) = {4, 1}; \n')
		
	# CREATION SURFACE DU RECTANGLE ##
	geo_file.write('Line Loop(5) = {10, 20, 30, 40} ;\n')
	geo_file.write('Plane Surface(6) = {5} ;\n')
		
	## MAILLAGE ##
	geo_file.write('Transfinite Line {10, 30} = Ceil(nx) Using Progression 1 ;\n')
	geo_file.write('Transfinite Line {20, 40} = Ceil(ny) Using Progression 1 ;\n')
	geo_file.write('Transfinite Surface {6}; \n')
	geo_file.write('Recombine Surface {6}; \n')
	geo_file.write('Coherence; \n')
		
	## DEFINITION LIGNE ET SURFACE PHYSIQUE ##
	geo_file.write('Physical Line(100) = {10}; \n')
	geo_file.write('Physical Line(200) = {20}; \n')
	geo_file.write('Physical Line(300) = {30}; \n')
	geo_file.write('Physical Line(400) = {40}; \n')
		
	geo_file.write('Physical Surface(1000) = {6}; \n')
	geo_file.close()
		
	os.system('gmsh '+self.geo_file+' -2 -order 1 -o '+self.geo_file[:-3]+'msh')
		
        with open('./MESH/'+self.name+'.msh', 'r') as file:
            data = file.readlines()
        
        ## NOMBRE DE NOEUDS ##
        ligne_nodes = 4
        self.nb_nodes = int(data[ligne_nodes])
        print "nb_nodes", self.nb_nodes
        self.ndof = 2*self.nb_nodes
 
        ## NOMBRE D'ELEMENTS ##
        ligne_elements = ligne_nodes+self.nb_nodes+3
        self.nelem = int(data[ligne_elements])
        
        ## COORDONNEES DES NOEUS ##        
        self.coor = np.empty((self.nb_nodes, 2))
        for k in range(ligne_nodes+1, ligne_nodes+1+self.nb_nodes):
			mots = string.split(data[k])
			self.coor[k-ligne_nodes-1,0] = float(mots[1])
			self.coor[k-ligne_nodes-1,1] = float(mots[2])

        ## LISTES ELEMENTS Q4 + NOEUDS ##
        self.nelm_q4 =0
        liste_elm_q4 = []
        for i in range(ligne_elements+1, ligne_elements+self.nelem+1):
			mots = string.split(data[i])
			if (mots[1] == "3"):
				liste_elm_q4.append(i)
				self.nelm_q4 = self.nelm_q4 +1
		
        self.elm_q4 = np.empty((self.nelm_q4, 4))
        k = 0
        for i in liste_elm_q4:
            mots = string.split(data[i])
            self.elm_q4[k, 0] = 	int(mots[-4])
            self.elm_q4[k, 1] = 	int(mots[-3])	
            self.elm_q4[k, 2] = 	int(mots[-2])
            self.elm_q4[k, 3] = 	int(mots[-1])
            k = k+1
					
  
            
    def get_nodes(self, group):
        nodes = []
		
        with open('./MESH/'+self.name+'.msh', 'r') as file:
            data = file.readlines()
        
        ## NOMBRE DE NOEUDS ##
        ligne_nodes = 4
        nb_nodes = int(data[ligne_nodes])
        
        ## NOMBRE D'ELEMENTS ##
        ligne_elements = ligne_nodes+nb_nodes+3
        nb_elements = int(data[ligne_elements])
        print 'nb_nodes : ', nb_nodes
        print 'nb_elm : ', nb_elements
        
        ## RECHERCHE DE NOEUD PAR GROUPE PHYSIQUE ##
        for i in range(ligne_elements+1, ligne_elements+nb_elements):
			mots = string.split(data[i])
			if (mots[3] == group)and(mots[1]=="1"):
				nodes.append(int(mots[-1])-1)
				nodes.append(int(mots[-2])-1)
        nodes = list(set(nodes))
        return nodes
        
            
    def form_q4(self, ksi, eta):
        "Affectation des fonctions de forme d'un element Q4"
        N = np.empty((2, 8))
		
        N1 = 0.25*(1. - ksi)*(1. - eta)
        N2 = 0.25*(1. + ksi)*(1. - eta)
        N3 = 0.25*(1. + ksi)*(1. + eta)
        N4 = 0.25*(1. - ksi)*(1. + eta)
		
        N[0,0] = N1
        N[0,1] = 0.
        N[0,2] = N2
        N[0,3] = 0.
        N[0,4] = N3
        N[0,5] = 0.
        N[0,6] = N4
        N[0,7] = 0.
		
        N[1,0] = 0.
        N[1,1] = N1
        N[1,2] = 0.
        N[1,3] = N2
        N[1,4] = 0.
        N[1,5] = N3
        N[1,6] = 0.
        N[1,7] = N4
		
        return N
		

    def dform_q4(self, ksi , eta):
		
        dN = np.empty((2, 4))
        
        dN1_ksi =   0.25*(-1.+eta)
        dN2_ksi =   0.25*(1.-eta)
        dN3_ksi =   0.25*(1.+eta)
        dN4_ksi =   0.25*(-1.-eta)
		
        dN1_eta =  0.25*(-1.+ksi)
        dN2_eta =  0.25*(-1.-ksi)
        dN3_eta =  0.25*(1.+ksi)
        dN4_eta =  0.25*(1.-ksi)
		
		
        dN[0, 0] = dN1_ksi
        dN[0, 1] = dN2_ksi
        dN[0, 2] = dN3_ksi
        dN[0, 3] = dN4_ksi

        dN[1, 0] = dN1_eta
        dN[1, 1] = dN2_eta
        dN[1, 2] = dN3_eta
        dN[1, 3] = dN4_eta
		
        #~ print 'dform_q4 :' +str(dN)
        return dN
		          
    def gaussPoint(self):
	"Definition des points de Gauss"
		
	self.Gauss_point_q4 = np.empty((4,2))
	p = 0.577350269189626
		
	self.Gauss_point_q4[0] = [-p, -p] 
	self.Gauss_point_q4[1] = [ p, -p] 
	self.Gauss_point_q4[2] = [ p,  p] 
	self.Gauss_point_q4[3] = [-p,  p] 
	
    def compute_D(self, E,v):
        self.D = np.empty((3,3))
        ## LIGNE 1
        self.D[0,0] = 1.
        self.D[0,1] = v
        self.D[0,2] = 0.
        ## LIGNE 2
        self.D[1,0] = v
        self.D[1,1] = 1.
        self.D[1,2] = 0.
        ## LIGNE 3
        self.D[2,0] = 0.
        self.D[2,1] = 0.
        self.D[2,2] = 0.5*(1.-v)
	
        self.D = self.D*(E/(1.-v**2.))
            
    def compute_D_plane_strain(self, E,v, J):
		
        u0 = E/(2*(1+v))
        l0 = v*E/((1+v)*(1-2.*v))
        
        l = l0
        u = u0 - l*np.log(J)
		
        self.D = np.empty((3,3))
        ## LIGNE 1
        self.D[0,0] = l+2.*u
        self.D[0,1] = l
        self.D[0,2] = 0.
        ## LIGNE 2
        self.D[1,0] = l
        self.D[1,1] = l+2.*u
        self.D[1,2] = 0.
        ## LIGNE 3
        self.D[2,0] = 0.
        self.D[2,1] = 0.
        self.D[2,2] = u
	
	
    def det(self,Matrix):
        a = np.linalg.det(Matrix)
        return a


    def compute_Bulk(self, GP,x, liste, nddl, emod, nue):
        Ki = np.zeros((8,8))
        K = np.zeros((nddl, nddl))
        t =1
        #~ print "noeuds : ", [liste[0]+1,liste[1]+1,liste[2]+1,liste[3]+1]

        ##LOCALISATION ##
        loc= [2*int(liste[0]), 2*int(liste[0])+1, 2*int(liste[1]), 2*int(liste[1])+1,2*int(liste[2]), 2*int(liste[2])+1,2*int(liste[3]), 2*int(liste[3])+1]

		## MATRICE DE RIGIDITE ELEMENTAIRE ##	
        for i in range(0, 4):
            N = self.form_q4(GP[i, 0], GP[i,1 ])
            dN = self.dform_q4(GP[i, 0], GP[i,1 ])     
            J = self.Jacob(dN, x)
            #~ print 'J : ', J
            det1 = self.det(J)
            #~ print 'det1 : ', det1
            self.compute_D_plane_strain(emod, nue, 1)
            B = self.B_matrix(GP[i,0],GP[i,1], J, dN)
            Ki = Ki + np.tensordot(np.tensordot(np.transpose(B), self.D, axes =1), B, axes =1)*det1
		
		## PRISE EN COMPTE DE L'EPAISSEUR t ##
        Ki = t*Ki
        
        ## ASSEMBLAGE MATRICE DE RIGIDITE GLOBALE ##
        for i in range(Ki.shape[0]):
			for j in range(Ki.shape[0]):
				K[int(loc[i]),int(loc[j])] = Ki[i, j]
        return K

    def Compute_global_stiffness_matrix(self, emod, nue):
	self.Initialize_stiffness()
	#~ self.compute_D(emod, nue)
	self.gaussPoint()
						## ITERATION SUR TOUT LES ELEMENTS ##
	for i in range(self.nelm_q4):
		liste_noeuds =  np.ndarray.tolist(self.elm_q4[i]-1)
		x = self.coor[liste_noeuds]
		K = self.compute_Bulk(self.Gauss_point_q4, x, liste_noeuds, self.ndof, emod, nue)
		self.K_tot = K + self.K_tot

    #~ def compute_fint(self, GP,x, liste,  P):
		
		#~ ## EPAISSEUR ELEMENTS ##
		#~ t =1
		
		#~ ## POINT DE GAUSS ##
		#~ GP = self.Gauss_point_q4
		
        #~ for i in range(self.nelm_q4):
			
			#~ ## NUMERO ET COOR NOEUDS ##
			#~ liste =  np.ndarray.tolist(self.elm_q4[i]-1)
			#~ x = self.coor[liste_noeuds]
			
			#~ ## INITIALISATION FINT ##
			#~ Fint = np.zeros((4,2))
			


			#~ ##LOCALISATION ##
			#~ loc= [2*int(liste[0]), 2*int(liste[0])+1, 2*int(liste[1]), 2*int(liste[1])+1,2*int(liste[2]), 2*int(liste[2])+1,2*int(liste[3]), 2*int(liste[3])+1]

			#~ ## MATRICE DE RIGIDITE ELEMENTAIRE ##	
			#~ for i in range(0, 4):
				#~ N = self.form_q4(GP[i, 0], GP[i,1 ])
				#~ dN = self.dform_q4(GP[i, 0], GP[i,1 ])     
				#~ J = self.Jacob(dN, x)
				#~ det1 = self.det(J)
				#~ B0 = self.B0(GP[i,0],GP[i,1], J, dN)
				#~ Fint = Fint + np.tensordot(np.transpose(B0), P, axes =1)*det1
				#~ print 'Fint :  ', Fint
		

    def Jacob(self,dN, x ):
		
        J = np.dot(dN, x)
        return J
		
    def B_matrix(self, ksi, eta, J, dN):
		
        B = np.empty((3,8))
		
        J_inv = np.linalg.inv(J)
		
		

        dN1_x = np.dot(J_inv, dN[:,0])[0]
        dN2_x = np.dot(J_inv, dN[:,1])[0]
        dN3_x = np.dot(J_inv, dN[:,2])[0]
        dN4_x = np.dot(J_inv, dN[:,3])[0]

        dN1_y = np.dot(J_inv, dN[:,0])[1]
        dN2_y = np.dot(J_inv, dN[:,1])[1]
        dN3_y = np.dot(J_inv, dN[:,2])[1]
        dN4_y = np.dot(J_inv, dN[:,3])[1]
		
        B[0] = [dN1_x, 0., dN2_x, 0. , dN3_x, 0., dN4_x, 0.]
        B[1] = [0., dN1_y, 0., dN2_y, 0., dN3_y, 0., dN4_y]
        B[2] = [dN1_y, dN1_x, dN2_y, dN2_x, dN3_y, dN3_x, dN4_y, dN4_x]
        
        return B

    def B_NON_LIN(self, ksi, eta, J, dN):
		
        B = np.empty((4,8))
		
        J_inv = np.linalg.inv(J)
		
		

        dN1_x = np.dot(J_inv, dN[:,0])[0]
        dN2_x = np.dot(J_inv, dN[:,1])[0]
        dN3_x = np.dot(J_inv, dN[:,2])[0]
        dN4_x = np.dot(J_inv, dN[:,3])[0]

        dN1_y = np.dot(J_inv, dN[:,0])[1]
        dN2_y = np.dot(J_inv, dN[:,1])[1]
        dN3_y = np.dot(J_inv, dN[:,2])[1]
        dN4_y = np.dot(J_inv, dN[:,3])[1]
		
        B[0] = [dN1_x, 0., dN2_x, 0. , dN3_x, 0., dN4_x, 0.]
        B[1] = [dN1_y, 0., dN2_y, 0. , dN3_y, 0., dN4_y, 0. ]
        B[2] = [0., dN1_x, 0., dN2_x, 0., dN3_x, 0., dN4_x]
        B[3] = [0., dN1_y, 0., dN2_y, 0., dN3_y, 0., dN4_y]
        
        return B

    def B0(self, ksi, eta, J, dN):
		
        B = np.empty((2,4))
		
        J_inv = np.linalg.inv(J)
		
		

        dN1_x = np.dot(J_inv, dN[:,0])[0]
        dN2_x = np.dot(J_inv, dN[:,1])[0]
        dN3_x = np.dot(J_inv, dN[:,2])[0]
        dN4_x = np.dot(J_inv, dN[:,3])[0]

        dN1_y = np.dot(J_inv, dN[:,0])[1]
        dN2_y = np.dot(J_inv, dN[:,1])[1]
        dN3_y = np.dot(J_inv, dN[:,2])[1]
        dN4_y = np.dot(J_inv, dN[:,3])[1]
		
        B[0] = [dN1_x, dN2_x, dN3_x, dN4_x]
        B[1] = [dN1_y, dN2_y, dN3_y, dN4_y]
        
        return B

    def plot_mesh(self):
        pl.figure()
        pl.plot(self.coor[:,0], self.coor[:,1], 'ob')
        pl.xlim(min(self.coor[:,0]) -10. , max(self.coor[:,0]) +10. )
        pl.ylim(min(self.coor[:,1]) -10. , max(self.coor[:,1]) +10. )
        pl.show()

    def plot_disp(self):
        pl.figure()
        pl.plot(self.coor[:,0], self.coor[:,1], 'ob')
        pl.plot(self.new_coor[:,0], self.new_coor[:,1], 'or')
        pl.xlim(min(self.new_coor[:,0]) -10. , max(self.new_coor[:,0]) +10. )
        pl.ylim(min(self.new_coor[:,1]) -10. , max(self.new_coor[:,1]) +10. )
        pl.show()

    def plot_U(self, k):
		X, Y = np.meshgrid(self.coor[:,0], self.coor[:,1])
		disp= scipy.interpolate.griddata((self.coor[:,0], self.coor[:,1]), self.U[:,k], (X, Y), method='linear')
		pl.figure()
		pl.pcolormesh(X, Y, disp,cmap='jet', vmin= min(self.U[:,k]), vmax=max(self.U[:,k]) )
		pl.colorbar()
		pl.show()
		
    def plot_E(self,E):
		X, Y = np.meshgrid(self.coor[:,0], self.coor[:,1])
		disp= scipy.interpolate.griddata((self.coor[:,0], self.coor[:,1]), E, (X, Y), method='linear')
		pl.figure()
		pl.pcolormesh(X, Y, disp,cmap='jet', vmin= min(E), vmax=max(E) )
		#~ pl.pcolormesh(X, Y, disp,cmap='jet' )
		pl.colorbar()
		pl.show()
        
    def Impose_displacement_by_penalisation(self, nodes, ddl , value, alpha):
	# FActeur terme penalisant
	alpha1 = 1e15*np.abs(np.max(echantillon.K_tot))
	for i in nodes:
		for j in ddl:
			self.F[2*int(i)+int(j)] = self.F[2*int(i)+int(j)] + alpha*value[j]
			self.K_tot[2*int(i)+int(j), 2*int(i)+int(j)] = self.K_tot[2*int(i)+int(j), 2*int(i)+int(j)] +alpha

    def Impose_displacement(self, nodes, ddl , value):
	for i in nodes:
		for j in ddl:
			#~ print 'yo'
			self.U_known_bool[2*int(i)+int(j)] = 1
			self.U_known[2*int(i)+int(j)] = value[j]
			self.F[2*int(i)+int(j)] = value[j]
			

    def Modify_Force_ext(self):
	
	LISTE_F = np.argwhere(self.U_known_bool ==0)
	LISTE_F = LISTE_F[:,0]

	LISTE_U = np.argwhere(self.U_known_bool ==1)
	LISTE_U = LISTE_U[:,0]
	
	for k in LISTE_F:
		for i in LISTE_U:
			if k != i:
				self.F[k] = self.F[k] - self.K_tot[k, i]*self.U_known[i]
				
    def Unit_diago_K(self):
	LISTE_U = np.argwhere(self.U_known_bool ==1)
	LISTE_U = LISTE_U[:,0]							
	for i in LISTE_U:
		self.K_tot[:, i] = 0.
		self.K_tot[i, :] = 0.
		self.K_tot[i, i] = 1.


    def Compute_localisation(self):
	Global_loc = []
	for i in range(self.nelm_q4):
		liste =  np.ndarray.tolist(self.elm_q4[i]-1)
		loc= [int(liste[0]), int(liste[1]), int(liste[2]),int(liste[3])]
		Global_loc.append(loc)
	return np.asarray(Global_loc)
	
    def Initialize_force(self):
	self.F = np.zeros((self.ndof, 1))

    def Initialize_stiffness(self):
	self.K_tot = np.zeros((self.ndof, self.ndof))

    def Compute_dof(self, u):
	self.new_coor = np.zeros((self.nb_nodes, 2))
	self.U = np.zeros((self.nb_nodes, 2))
	for i in range(self.coor.shape[0]):
		self.new_coor[i, 0] = self.coor[i, 0] + u [2*i]
		self.new_coor[i, 1] = self.coor[i, 1] + u [2*i +1]
		
		self.U[i, 0] = u[2*i]
		self.U[i, 1] = u[2*i+1]

    def Compute_GL_deformations(self, emod, nue):
	E11 = []
	E12 = []
	E21 = []
	E22 = []
	
	S11 = []
	S22 = []
	S21 = []
	
	Fint_test = np.zeros((self.ndof, 1))
						## ITERATION SUR TOUT LES ELEMENTS ##
	for i in range(self.nelm_q4):
		liste_noeuds =  np.ndarray.tolist(self.elm_q4[i]-1)
		x = self.coor[liste_noeuds]
		Ei11, Ei12,Ei21,Ei22, Fint, loc = self.Compute_gradient_and_stress(self.Gauss_point_q4, x, liste_noeuds, self.U[liste_noeuds, :], emod, nue)
		
                Fint_test[loc[0]] = Fint_test[loc[0]] +  Fint[0, 0]
                Fint_test[loc[1]] = Fint_test[loc[1]] +  Fint[0, 1]
                Fint_test[loc[2]] = Fint_test[loc[2]] +  Fint[1, 0]
                Fint_test[loc[3]] = Fint_test[loc[3]] +  Fint[1, 1]
                Fint_test[loc[4]] = Fint_test[loc[4]] +  Fint[2, 0]
                Fint_test[loc[5]] = Fint_test[loc[5]] +  Fint[2, 1]
                Fint_test[loc[6]] = Fint_test[loc[6]] +  Fint[3, 0]
                Fint_test[loc[7]] = Fint_test[loc[7]] +  Fint[3, 1]
                
        
		#~ stop_here
		E11.append(Ei11)
		E12.append(Ei12)
		E21.append(Ei21)
		E22.append(Ei22)
	print Fint_test
	return E11, E12, E21, E22, Fint_test
	
    def Compute_gradient_and_stress(self, GP, x, liste_noeuds, U, emod, nue):
	Ei11 = []
	Ei12 = []
	Ei21 = []
	Ei22 = []
	P_glob = []
	## INITIALISATION FINT ##
	Fint = np.zeros((4,2))
	
	for i in range(0, 4):
		
			## DEFORMATION GRADIENT ##
			## VOIR : http://www.continuummechanics.org/finiteelementmapping.html ##
            dN = self.dform_q4(GP[i, 0], GP[i,1 ]) 
            #~ print 'dN : ', dN
            #~ print 'U  : ', U
            Fi = np.dot(np.dot(dN, U), np.linalg.inv(np.dot(dN,x))) + np.eye(2,2)
            #~ print 'Fi : ', Fi

            
            ## GL STRAINS ##
            ## POSSIBILITE D'OBTENIR E SANS PASSER PAR F : VOIR FEM, BATHE p 570 ##
            Ei = 0.5*(np.dot(Fi.T, Fi)-np.eye(2,2))
            Ei11.append(Ei[0,0])
            Ei12.append(Ei[0,1])
            Ei21.append(Ei[1,0])
            Ei22.append(Ei[1,1])
            
            ## SECOND PIOLA KIRCHOFF STRESS ##
            J = np.linalg.det(Fi)
            self.compute_D_plane_strain(emod, nue, J)  
            S = np.dot(self.D, [Ei[0,0], Ei[1,1], 2.*Ei[1,0]] )
            Si = np.zeros((2, 2))
            Si[0, 0] = S[0]
            Si[0, 1] = S[2]
            Si[1, 0] = S[2]
            Si[1, 1] = S[1]
            
            ## FIRST PIOLA KIRCHOFF STRESS ##
            P = np.dot(Si, Fi.T)

					## FORCE INTERIEURE ##
            ## NUMERO ET COOR NOEUDS ##
            #~ liste =  np.ndarray.tolist(self.elm_q4[i]-1)
            x = self.coor[liste_noeuds]
            
			##LOCALISATION ##
            loc= [2*int(liste_noeuds[0]), 2*int(liste_noeuds[0])+1, 2*int(liste_noeuds[1]),2*int(liste_noeuds[1]) +1 ,2*int(liste_noeuds[2]), 2*int(liste_noeuds[2])+1,2*int(liste_noeuds[3]), 2*int(liste_noeuds[3])+1]
			
            J = self.Jacob(dN, x)
            det1 = self.det(J)
            B0 = self.B0(GP[i,0],GP[i,1], J, dN)
            Fint = Fint + np.tensordot(np.transpose(B0), P, axes =1)*det1
				            

        return Ei11, Ei12, Ei21, Ei22, Fint, loc

    def Extrapolate_at_node(self, GP, E11, E12, E22):
	## TRANSFORMATION DES PGs pour EXTRAPOLATION ##
	GP = 3.*GP
	
	## TRANSFORMATION DES LISTES EN TABLEAUX ##
	E11 = np.asarray(E11)
	E12 = np.asarray(E12)
	E22 = np.asarray(E22)
	
	## CREATION DE NOUVEAUX TABLEAUX ##
	E11_nodes = np.zeros((E11.shape[0], E11.shape[1]))
	E21_nodes = np.zeros((E12.shape[0], E12.shape[1]))
	E22_nodes = np.zeros((E22.shape[0], E22.shape[1]))
	
	
	EXTRA = np.zeros((4 , 4))
	for i in range(GP.shape[0]):
		EXTRA[i, 0] = 0.25*(1. - GP[i,0])*(1. - GP[i,1])
		EXTRA[i, 1] = 0.25*(1. + GP[i,0])*(1. - GP[i,1])
		EXTRA[i, 2] = 0.25*(1. + GP[i,0])*(1. + GP[i,1])
		EXTRA[i, 3] = 0.25*(1. - GP[i,0])*(1. + GP[i,1])
	
	for i in range(E11.shape[0]):
		E11_nodes[i] = np.dot(EXTRA, E11[i])
		E21_nodes[i] = np.dot(EXTRA, E12[i])
		E22_nodes[i] = np.dot(EXTRA, E22[i])
	return E11_nodes, E22_nodes, E21_nodes
	
		
    def Averaging_at_node(self, E11_nodes, E22_nodes, E12_nodes, LOC):
        ## MOYENNE DES DEFORMATIONS AUX NOEUDS ##
        E11_mean = np.zeros((self.nb_nodes,1))
        E12_mean = np.zeros((self.nb_nodes,1))
        E22_mean = np.zeros((self.nb_nodes,1))

        for i in range(self.nb_nodes):
            loc = np.argwhere(LOC == i)
            for k in range(loc.shape[0]):
                E11_mean[i] = E11_mean[i] + E11_nodes[loc[k,0], loc[k,1]]
                E22_mean[i] = E22_mean[i] + E22_nodes[loc[k,0], loc[k,1]]
                E12_mean[i] = E12_mean[i] + E12_nodes[loc[k,0], loc[k,1]]
            E11_mean[i] = E11_mean[i]/float(loc.shape[0])
            E22_mean[i] = E22_mean[i]/float(loc.shape[0])
            E12_mean[i] = E12_mean[i]/float(loc.shape[0])
        return E11_mean, E22_mean,E12_mean
