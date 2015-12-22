#!/usr/bin/python

import os, numpy, glob, math, collections, linecache, itertools, tkMessageBox
import matplotlib.pyplot as plt
from tkFileDialog   import askopenfilename
from Tkinter import *

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#Nf = 60
Lem = 0.04125
A = 0.8*0.042

filename = "tirant_macro_CN2"
filename0 = filename + ".data"
filename1 = filename + ".gid.msh"
filename2 = filename + ".gid.res"
filename3 = filename + ".reac"

for J in range(1,5):

	os.chdir("/home/christian/Documents/FIDES/FIDES_POST/FIDES_resultats/" + filename + "/" + str(J))

	f = numpy.loadtxt(filename3)
	lenmax = 0
	for i in f:
		if str(i[0]) == "nan":
			break
		else:
			lenmax += 1

	pas_min = 0
	pas_max = int(lenmax)
	pas_int = 1

	if pas_max != 301:
		continue

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
	#surf_mid = 0.5*(surf_sup + surf_inf)

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
    
	#ln = min([i[0] for i in noeuds_mid.values() if i[0] != x_mid_min]) - x_mid_min
	ln = min([i[0] for i in noeuds_inf.values() if i[0] != x_inf_min]) - x_inf_min
	#Le = ln * round(L_mid/Nf/ln)
	Le = ln * round(Lem/ln)
	Le = Lem
	#iLe = [i*Le + x_mid_min for i in range(int(round(L_mid/Le))+1)]
	iLe = [i*Le + x_inf_min for i in range(int(round(L_inf/Le))+1)]

	faces = {}
	for j in iLe:
		for i in [k for k in coor if abs(coor[k][0]-j) < 0.00001]:
			if coor[i][0] not in faces.keys():
				faces[coor[i][0]] = []
			else:
				faces[coor[i][0]].append(i)

	faces = collections.OrderedDict(sorted(faces.items()))

	noeuds_dep = {}
	force = []
	for z in range(pas_min, pas_max, pas_int):
		o = linecache.getline(filename3, z +2)
		force.append(float("".join([o[k] for k in range(75,104)])))
		start_dep = 3 + 3*z + 4*z + 2*z*noeuds
		for i in range(noeuds):
			m = linecache.getline(filename2, start_dep + i + 1)
			j = int("".join([m[k] for k in range(12)]))
			if j in list(itertools.chain.from_iterable(faces.values())):
				if j not in noeuds_dep.keys():
					noeuds_dep[j] = []
				noeuds_dep[j].append(float("".join([m[k] for k in range(13,38)])))

	stress = [i/A for i in force]

	dep = {}
	delta = {}
	eps = {}

	for i in range(len(faces)):
		dep[i] = []
		for k in faces[faces.keys()[i]]:
			dep[i].append(noeuds_dep[k])
		valsave = [numpy.mean(list(zip(*dep[i])[j])) for j in range(len(dep[i][0]))]
		dep[i] = valsave
		if i>=1:
			if abs(faces.keys()[i]-faces.keys()[i-1] - Le) < 0.00001:
				delta[i] = [dep[i][j] - dep[i-1][j] for j in range(len(dep[0]))]
				eps[i] = [j/Le for j in delta[i]]
			else:
				continue

	linecache.clearcache()

	os.chdir("/home/christian/Desktop/coupes_" + filename + "/force-dep")

	for i in delta:
		f = open(filename0[:-5] + "_force-depl_t-" + str(J) + "_e-" + str(i) + ".txt", "w")
		for j in range(len(delta[i])):
			f.write('%.17E' % delta[i][j] + "  ")
			f.write('%.17E' % force[j] + "\n")
		f.close()

	os.chdir("/home/christian/Desktop/coupes_" + filename + "/sig-eps")

	for i in eps:
		f = open(filename0[:-5] + "_sig-eps_t-" + str(J) + "_e-" + str(i) + ".txt", "w")
		for j in range(len(eps[i])):
			f.write('%.17E' % eps[i][j] + "  ")
			f.write('%.17E' % stress[j] + "\n")
		f.close()
    
os.chdir("/home/christian/Documents/MyThesis/PyFun")
