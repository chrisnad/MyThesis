#!/usr/bin/python

import os, numpy, glob, math, collections, linecache, itertools
import matplotlib.pyplot as plt

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#Nf = 60
Lem = 0.04125
A = 0.8*0.042

pas_min = 0
pas_max = 301
pas_int = 1

os.chdir("/home/christian/Documents/FIDES/FIDES_POST_test/FIDES_resultats/tirant_macro_CN/1")

filename = "tirant_macro_CN"
filename0 = filename + ".data"
filename1 = filename + ".gid.msh"
filename2 = filename + ".gid.res"
filename3 = filename + ".reac"

noeuds = int(''.join([k for k in linecache.getline(filename0, 16) if k.isdigit()]))

coor = {}
for i in range(noeuds):
    a = ""
    for j in linecache.getline(filename1, i+3):
        if j.isdigit() or j == "." or j == "E" or j == "-":
            a += j
        else:
            if a.isdigit():
                b = int(a)
                coor[b]=[]
            elif is_number(a):
                coor[b].append(float(a))
            a = ""

surf_sup = max([i[1] for i in coor.values()])
surf_inf = min([i[1] for i in coor.values()])
surf_mid = 0.5*(surf_sup + surf_inf)

noeuds_sup = {}
noeuds_inf = {}
#noeuds_mid = {}

for i in [j for j in coor if coor[j][1] == surf_sup]:
    noeuds_sup[i] = coor[i]
for i in [j for j in coor if coor[j][1] == surf_inf]:
    noeuds_inf[i] = coor[i]
#for i in [j for j in coor if coor[j][1] == surf_mid]:
    #noeuds_mid[i] = coor[i]

x_sup_max = max([i[0] for i in noeuds_sup.values()])
x_sup_min = min([i[0] for i in noeuds_sup.values()])
x_inf_max = max([i[0] for i in noeuds_inf.values()])
x_inf_min = min([i[0] for i in noeuds_inf.values()])
#x_mid_max = max([i[0] for i in noeuds_mid.values()])
#x_mid_min = min([i[0] for i in noeuds_mid.values()])

L_sup = x_sup_max - x_sup_min
L_inf = x_inf_max - x_inf_min
#L_mid = x_mid_max - x_mid_min

ln = min([i[0] for i in noeuds_inf.values() if i[0] != x_inf_min]) - x_inf_min
#Le = ln * round(L_mid/Nf/ln)
Le = ln * round(Lem/ln)
Le = Lem
iLe = [i*Le + x_inf_min for i in range(int(round(L_inf/Le))+1)]

faces = {}
for j in iLe:
    for i in [k for k in coor if abs(coor[k][0]-j) < 0.00001]:
        if coor[i][0] not in faces.keys():
            faces[coor[i][0]] = []
        else:
            faces[coor[i][0]].append(i)

faces = collections.OrderedDict(sorted(faces.items()))

#noeuds_cont = {}
noeuds_dep = {}
force = []
for z in range(pas_min, pas_max, pas_int):
    o = linecache.getline(filename3, z +2)
    force.append(float("".join([o[k] for k in range(75,104)])))
    start_dep = 3 + 3*z + 4*z + 2*z*noeuds
    #start_cont = 7 + noeuds + 3*z + 4*z + 2*z*noeuds
    for i in range(noeuds):
        m = linecache.getline(filename2, start_dep + i + 1)
        j = int("".join([m[k] for k in range(12)]))
        if j in list(itertools.chain.from_iterable(faces.values())):
            if j not in noeuds_dep.keys():
                noeuds_dep[j] = []
            noeuds_dep[j].append(float("".join([m[k] for k in range(13,38)])))
        #n = linecache.getline(filename2, start_cont + i + 1)
        #j = int("".join([n[k] for k in range(12)]))
        #if j in list(itertools.chain.from_iterable(faces.values())):
            #if j not in noeuds_cont.keys():
                #noeuds_cont[j] = []
            #noeuds_cont[j].append(float("".join([n[k] for k in range(13,38)])))

stress = [i/A for i in force]

dep = {}
delta = {}
eps = {}
#cont = {}

for i in range(len(faces)):
    dep[i] = []
    #cont[i] = []
    for k in faces[faces.keys()[i]]:
        dep[i].append(noeuds_dep[k])
        #cont[i].append(noeuds_cont[k])
    dep[i] = [numpy.mean(list(zip(*dep[i])[j])) for j in range(len(dep[i][0]))]
    #cont[i] = [numpy.mean(list(zip(*cont[i])[j])) for j in range(len(cont[i][0]))]
    if i>=1:
		if abs(faces.keys()[i]-faces.keys()[i-1] - Le) < 0.00001:
			delta[i] = [dep[i][j] - dep[i-1][j] for j in range(len(dep[0]))]
			eps[i] = [j/Le for j in delta[i]]
		else:
			continue
			
linecache.clearcache()

f1 = plt.figure(1)
ax1 = f1.add_subplot(111)
for i in range(1 ,len(delta)+1):
    ax1. plot(delta[i], force)
plt.title("Force / Depl")

f2 = plt.figure(2)
ax2 = f2.add_subplot(111)
for i in range(1, len(eps)+1):
    ax2.plot(eps[i], stress)
plt.title("Sigma / Eps")

#f3 = plt.figure(3)
#ax3 = f3.add_subplot(111)
#for i in range(1, len(eps)+1):
#    ax3.plot(eps[i], cont[i])
#plt.title("Sigma / Eps")

plt.show()
    
