# -*-coding: utf-8 -*-
from nlin_fem import *
#~ import nlin_fem as nlin
import pylab as pl
import numpy as np
import scipy.interpolate





###########################################################################

					## MISE EN DONNEES DU PrOBLEME ##

###########################################################################
			## PREPARATION DE L'ECHANTILLON / MAILLAGE  ##			
echantillon =sample(nom = "yo",materiau = "mou",maillage = "quad")
echantillon.mesh(lx = 10., ly = 10., nx = 15., ny=15.)
echantillon.plot_mesh()

					## TABLE CONNECTIVITE TOTALE ##
LOC = echantillon.Compute_localisation()

			## EXTRACTION DES NOEUDS PAR GROUPE PHYSIQUE POUR CL##			
## LIGNE  INF ##
nodes_inf = echantillon.get_nodes(group = "100")
## LIGNE SUP ## 
nodes_sup = echantillon.get_nodes(group = "300")
## LIGNE DROITE ##
nodes_d = echantillon.get_nodes(group = "200")


###########################################################################
						## PARAMETRES MATERIAU ##
					# MODELE NEO-HOOKEAN COMPRESSIBLE #
	#(J.BONNET, equation 6.27 : Ψ=0.5*μ(I C − 3) − μ ln J + 0.5*λ(ln J) 2) #
					 #S = μ(I − Cinv ) + λ(ln J)Cinv #
					 #σ =(μ/J)(b − I) + (λ/J)(ln J)I #					 
emod = 3. # YOUng
nue = 0.45 # Coeff poisson
						## PARAMETRES MATERIAU ##
###########################################################################


###########################################################################

						## SOFT NON-LINEAR MECHANICS ##

###########################################################################
							## MATRICE FORCE ##
echantillon.Initialize_force()

							## CONDITION AUX LIMITES ##
echantillon.U_known_bool = np.zeros((echantillon.ndof, 1))
echantillon.U_known= np.zeros((echantillon.ndof, 1))


				  ## APPLICATION DES CONDITIONS AUX LIMITES ##	
echantillon.Impose_displacement(nodes = nodes_inf, ddl = [0,1], value = [0.,0.] )
echantillon.Impose_displacement(nodes = nodes_sup, ddl = [0,1], value = [0,0.1])




							## STIFNESS MATRIX ##
echantillon.Compute_global_stiffness_matrix(emod, nue)

echantillon.Modify_Force_ext()
echantillon.Unit_diago_K()

						## CALCUL DES DEPLACEMENT ##
u = np.linalg.solve(echantillon.K_tot, echantillon.F)
echantillon.Compute_dof(u)


						## PLOT DES DEPLACEMENTS ##
echantillon.plot_disp()
## Disp / x ##
echantillon.plot_U(0)
## Disp / y ##
echantillon.plot_U(1)


						## CALCUL DES DEFORMATIONS ##
## DEFORMATION AUX PGs ##
E11, E12, E21, E22, Fint = echantillon.Compute_GL_deformations(emod, nue)
print 'E11[0] :  ', E11[0]
						## EXTRAPOLATION AU NOEUDS ##
E11_nodes, E22_nodes, E12_nodes = echantillon.Extrapolate_at_node(echantillon.Gauss_point_q4,E11, E12, E22 )

E11_mean, E22_mean, E12_mean = echantillon.Averaging_at_node(E11_nodes, E22_nodes, E12_nodes , LOC)
						## PLOT DES DEFORMATIONS ##
echantillon.plot_E(E11_mean)
echantillon.plot_E(E22_mean)
echantillon.plot_E(E12_mean)


print 'Fint-Fext : ', Fint-echantillon.F
## CALCUL DU RESIDU ##
# Force interieure
# Matrice tangente





