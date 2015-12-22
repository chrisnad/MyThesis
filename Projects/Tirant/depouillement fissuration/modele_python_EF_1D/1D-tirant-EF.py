import math, os, random
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
import sys

###     Géométrie
#       Longueur du tirant (m)

L = 1.75

#       Densité de fissures (/m)

den_fiss = 30/L

#       Section (m^2)

h = 0.072
ep = 0.8
A = h*ep

#       Section d'acier de diametre D (n_s aciers)

D_s = 0.012
n_s = 5
A_s = n_s*0.25*math.pi*D_s**2

###
#       Chargement appliqué : Déplacement imposé en x=L (en m)

Uimp = 0.00005
UimpEvol = [0.00005*i*0.1 for i in range(10) ] + [0.00005*(i+1) for i in range(200)]

###    Discrétisation de la géométrie
#       Nombre d'éléments sur la longueur

nelt = 50
nno = nelt + 1                                  #       Nombre de noeuds du maillage
Le = L/nelt                                       #       Longueur d'un élément
xe = [Le*i for i in range(nelt+1)]          #       Coordonnées des noeuds du maillage

#       Résistance à la traction des éléments

dataExp1 = open('n_fiss.txt', 'r')

L = [i for i in dataExp1]
Xf = []
Yf = []
for j in L:
        a=""
        b=""
        for i in range(22):
                a+=j[i]
                b+=j[i+25]
        Xf.append(float(a))
        Yf.append(float(b))

v = list(np.random.uniform(0.0, 1.0, nelt))

Xf.pop(0)
Yf.pop(0)

xp = Xf
yp = [i/max(Yf) for i in Yf]

f = interp1d(yp, xp, kind='cubic')
F = list(f(v))

RT = [i/A for i in F]
RTmin = min(RT)
RTmax = max(RT)
alp = [(i-RTmin)/(RTmax-RTmin) for i in RT]

###
#       Probabilité d'avoir un élément qui fissure

p_fiss = den_fiss * Le
rupt_el = list(np.random.binomial(1, p_fiss, nelt))
rt = []
for i in range(len(RT)):
        rt.append((1-rupt_el[i])*1e10 + rupt_el[i]*RT[i])
RT = rt
n_fiss = sum(rupt_el)
pn_fiss = float(n_fiss)/nelt

###                                     Materiaux
#        Acier (MPa)
E_s = 200000
#        Beton (MPa)
E_b = 35000
#        Equivalent A/B
E_m = (A_s*E_s+(A-A_s)*E_b)/A
#        Equivalent B_a fissure         (vecteur)
E_ms = (A_s*E_s)/A
E_mfmoy = p_fiss*E_m*E_ms/(E_m-(1-p_fiss)*E_ms)
E_min = pn_fiss*E_m*E_ms/(E_m-(1-pn_fiss)*E_ms)
E_mf = [(1.0/(1.0/E_min + alp[i] * (1.0/E_mfmoy - 1.0/E_min)))*rupt_el[i] for i in range(len(alp))]


###     MEF

ndlt = 1*nno
U = [0 for i in range(ndlt)]
resu = {'depl': [], 'cont': [], 'fiss': []}
FX = []
FY = []

#       Boucle sur nombre de pas de calculs

npas = len(UimpEvol)
lam = [1.0 for i in range(npas)]
Dini = [0.0 for i in range(nelt)]
D = Dini

### checked up to here ###

for ipas in range(npas):

        #       Initialisation pour chaque pas
        
        nconv = False
        dUt = [0.0 for i in range(ndlt)]
        print (" Pas de calcul n°: " + str(ipas) + "    Coefficient de charge " + str(lam[ipas]))

        #       Boucle sur le nombre de pas de calculs       
               
        for itera in range(5000):

                #       Initialisation du système matriciel
               
                K = [[0.0 for i in range(ndlt)] for j in range(ndlt)]
                R = [0.0 for i in range(ndlt)]
                Fint = [0.0 for i in range(ndlt)]
                sig = [0.0 for i in range(nelt)]

                #       boucle sur le nbr d'elements

                for ie in range(nelt):

                        #       Calcul et integration des contraintes sur l'element

                        eps = (U[ie+1]-U[ie])/Le
                        S = A
                        E = E_m

                        D[ie] = Dini[ie]
                        if eps > RT[ie]/E_m:
                                D[ie] = 1 - RT[ie]/(E_m*eps)
                                if eps > RT[ie]/E_mf[ie]:
                                        D[ie] = 1 - E_mf[ie]/E_m
                                
                                D[ie] = max([Dini[ie], D[ie]])
                                D[ie] = min(1.0, D[ie])

                        sig[ie] = (1 - D[ie])*E*eps
                        Fint_e = [sig[ie]*S*i for i in [-1, 1]]
                        Fint[ie] , Fint[ie+1] = Fint[ie] + Fint_e[0], Fint[ie+1] + Fint_e[1]

                        #       Assemblage du vec Residu

                        R[ie], R[ie+1] = R[ie] + Fint_e[0], R[ie+1] + Fint_e[1]

                        #       Assemblage matrice de rigidite

                        K[ie][ie], K[ie+1][ie] = K[ie][ie] + (1-D[ie])*E*S/Le, K[ie+1][ie] - (1-D[ie])*E*S/Le
                        K[ie][ie+1], K[ie+1][ie+1] = K[ie][ie+1] - (1-D[ie])*E*S/Le, K[ie+1][ie+1] + (1-D[ie])*E*S/Le

                #       CL

                if itera == 0 :
                        r = [R[i] - K[:][ndlt-1][i]*lam[ipas]*(-Uimp) for i in range(len(R))]
                        R = r
                        R[ndlt-1] = -lam[ipas]*Uimp
                else:
                        R[ndlt-1] = 0

                K[0] = [0.0 for i in range(ndlt)]
                for i in range(ndlt):
                        K[i][0] = 0.0
                K[0][0] = 1.0
                R[0] = 0.0

                K[ndlt-1] = [0.0 for i in range(ndlt)]
                for i in range(ndlt):
                        K[i][ndlt-1] = 0.0
                K[ndlt-1][ndlt-1] = 1.0

                #       Resolution du systeme matriciel

                dU = list(np.dot(np.linalg.inv(-np.array(K)), R))

                dUt = [dUt[i] + dU[i] for i in range(len(dU))]
                U = [U[i] + dU[i] for i in range(len(U))]

                #       test de convergence

                ndu = np.linalg.norm(dU)/(np.linalg.norm(dUt) + 1e-20)
                print (" Itération " + str(itera) + "    -       Critere : " + str(ndu))
                if ndu < 1e-4:
                        nconv = True
                        print ("        ")
                        break

        if nconv == False:
                break

        #       Stockage des resultats a chaque pas

        resu["depl"].append(U)
        resu["cont"].append(sig)
        resu["fiss"].append(sum(x > 0 for x in D))

        FX.append(U[ndlt-1])
        FY.append(Fint[ndlt-1])
                                   
f1 = plt.figure(1)
ax1 = f1.add_subplot(111)
ax1.plot(FX, FY)
plt.show()                        


