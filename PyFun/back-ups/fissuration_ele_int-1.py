import os, numpy, pylab, glob, math, matplotlib #sudo apt-get install python-matplotlib
from Tkinter import *
from tkFileDialog   import askopenfilename 
import tkMessageBox
import itertools

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

directory0 = "/home/cnader/Documents/FIDES/FIDES_POST/FIDES_resultats/test-list"
filename0 = "test-list.data"
directory1 = "/home/cnader/Documents/FIDES/FIDES_POST/FIDES_resultats/test-list"
filename1 = "test-list.gid.msh"
directory2 = "/home/cnader/Documents/FIDES/FIDES_POST/FIDES_resultats/test-list"
filename2 = "test-list.gid.res"

pas_min = 4
pas_max = 11
pas_int = 2

os.chdir(directory0)
f = open(filename0, 'r')

count = 0
for i in f:
    count += 1
    if count == 16:
        a = ""
        for j in i:
            if j.isdigit():
                a += j
            else:
                if a.isdigit():
                    noeuds = int(a)
                    break
        break

os.chdir(directory1)
f = open(filename1, 'r')

count = 0
coor = {}
for i in f:
    count += 1
    if count > 2 and count < 3 + noeuds:
            a = ''
            for j in i:
                if j.isdigit() or j == "." or j == "E" or j == "-":
                    a += j
                else:
                    if a.isdigit():
                        b = int(a)
                        coor[b]=[]
                    elif is_number(a):
                        coor[b].append(float(a))
                    a = ''
                    
surf_sup = -float("inf")
surf_inf = float("inf")

for i in coor:
    if coor[i][1] > surf_sup:
        surf_sup = coor[i][1]
    if coor[i][1] < surf_inf:
        surf_inf = coor[i][1]
       
noeuds_sup = {}
noeuds_inf = {}

for i in coor:
	if coor[i][1] == surf_sup:
		noeuds_sup[i] = coor[i]
	elif coor[i][1] == surf_inf:
		noeuds_inf[i] = coor[i]
		
noeuds_sup_fiss = []
noeuds_inf_fiss = []

while len(noeuds_sup) > 0:
	i = noeuds_sup.keys()[0]
	a = []
	for k in noeuds_sup:
		if k != i:
			if noeuds_sup[i] == noeuds_sup[k]:
				a.append(i)
				a.append(k)
	a = list(set(a))
	for j in a:
		del noeuds_sup[j]
	if len(a) != 0:
		noeuds_sup_fiss.append(a)
	else:
		del noeuds_sup[i]

while len(noeuds_inf) > 0:
	i = noeuds_inf.keys()[0]
	a = []
	for k in noeuds_inf:
		if k != i:
			if noeuds_inf[i] == noeuds_inf[k]:
				a.append(i)
				a.append(k)
	a = list(set(a))
	for j in a:
		del noeuds_inf[j]
	if len(a) != 0:
		noeuds_inf_fiss.append(a)
	else:
		del noeuds_inf[i]

os.chdir(directory2)

courbe_nbr_inf = []
courbe_nbr_sup = []
courbe_ouv_inf = []
courbe_ouv_sup = []
courbe_esp_inf = []
courbe_esp_sup = []

for z in range(pas_min, pas_max, pas_int):
    count = 0
    disp = []
    start = 3 + 3*z + 4*z + 2*z*noeuds              #equation change if it's 3D
    finish = 3 + 3*z + 4*z + (2*z + 1)*noeuds       #equation change if it's 3D
    f = open(filename2, 'r')
    for i in f:
        count += 1
        if count > start and count <= finish:
            for j in noeuds_inf_fiss:
                if (count - start) in j:
                    disp.append(i)
                    break
            for j in noeuds_sup_fiss:
                if (count - start) in j:
                    disp.append(i)
                    break
    
    noeuds_disp = {}
    for i in disp:
        a = ''
        for j in i:
            if j.isdigit() or j == "." or j == "E" or j == "-":
                a += j
            else:
                if a.isdigit():
                    b = int(a)
                    noeuds_disp[b]=[]
                elif is_number(a):
                    noeuds_disp[b].append(float(a))
                a = ''

    ouv_fiss_sup = []
    ouv_fiss_inf = []
    coo_fiss_sup = []
    coo_fiss_inf = []
    nbr_fiss_sup = 0
    nbr_fiss_inf = 0
    seuil = 3.0e-5

    for i in noeuds_inf_fiss:
        fiss = 0
        maxi = 0
        j = list(itertools.combinations(i, 2))
        for k in j:
            fiss = math.sqrt((noeuds_disp[k[0]][0] - noeuds_disp[k[1]][0])**2 + (noeuds_disp[k[0]][1] - noeuds_disp[k[1]][1])**2)
            if fiss > maxi:
                maxi = fiss
        if maxi >= seuil:
            ouv_fiss_inf.append(maxi)
            nbr_fiss_inf += 1
            coo_fiss_inf.append(coor[i[0]][0])

    for i in noeuds_sup_fiss:
        fiss = 0
        maxi = 0
        j = list(itertools.combinations(i, 2))
        for k in j:
            fiss = math.sqrt((noeuds_disp[k[0]][0] - noeuds_disp[k[1]][0])**2 + (noeuds_disp[k[0]][1] - noeuds_disp[k[1]][1])**2)
            if fiss > maxi:
                maxi = fiss
        if maxi >= seuil:
            ouv_fiss_sup.append(maxi)
            nbr_fiss_sup += 1
            coo_fiss_sup.append(coor[i[0]][0])

    if len(ouv_fiss_sup) > 0:
        moy_fiss_sup = sum(ouv_fiss_sup) / float(len(ouv_fiss_sup))
    else:
        moy_fiss_sup = 0

    if len(ouv_fiss_inf) > 0:
        moy_fiss_inf = sum(ouv_fiss_inf) / float(len(ouv_fiss_inf))
    else:
        moy_fiss_inf = 0
        
    courbe_nbr_inf.append(nbr_fiss_inf)
    courbe_nbr_sup.append(nbr_fiss_sup)
    courbe_ouv_inf.append(moy_fiss_inf)
    courbe_ouv_sup.append(moy_fiss_sup)
    
