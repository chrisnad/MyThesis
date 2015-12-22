#!/usr/bin/python

import numpy, os, glob, pylab
from scipy.optimize import fmin
import matplotlib.pyplot as plt
from shapely.geometry.polygon import  LinearRing, Polygon
from collections import *
#from statistics import *

filename1 = "tirant_fin_train_pl"
wd1 = "/home/christian/Documents/FIDES/FIDES_Resultats_calculs/Results/" + filename1
filename2 = "tirant_macro_CN2" #model
wd2 = "/home/christian/Documents/FIDES/FIDES_Resultats_calculs/Results/" + filename2

colonne_1 = 1
colonne_2 = 4

##################################

def area(p):
	return 0.5*abs(sum(x0*y1 - x1*y0 for ((x0, y0), (x1,y1)) in segments(p)))

def segments(p):
	return zip(p, p[1:]+[p[0]])

os.chdir(wd1)
Files = filename1 + "*"

mean_Y = {}

for fichier in glob.glob(Files):
	datalist = numpy.loadtxt(fichier)
	n = len(datalist[0])
	data = []
	for i in range(n):
		data.append(datalist[:,i])
	X = [x for x in data[colonne_1 - 1] if str(x)!="nan"]
	Y = [y for y in data[colonne_2 - 1] if str(y)!="nan"]
	if colonne_1 == "4":
		XX = [ float('%.2f' % elem) for elem in X ]
	else:
		XX = [ float('%.8f' % elem) for elem in X ]
	for i in range(len(XX)):
		if XX[i] in mean_Y:
			mean_Y[XX[i]].append(Y[i])
		else:
			mean_Y[XX[i]] = []
			mean_Y[XX[i]].append(Y[i])
			
for i in mean_Y:
	mean_Y[i] = sum(mean_Y[i])/len(mean_Y[i])

mean_Y1 = OrderedDict(sorted(mean_Y.items(), key=lambda t: t[0]))

##################################

os.chdir(wd2)
Files = filename2 + "*"

xy1 = list(zip(mean_Y1.keys(), mean_Y1.values()))
f = []

for fichier in glob.glob(Files):
	datalist = numpy.loadtxt(fichier)
	n = len(datalist[0])
	data = []
	dict_Y = {}
	for i in range(n):
		data.append(datalist[:,i])
	X = [x for x in data[colonne_1 - 1] if str(x)!="nan"]
	Y = [y for y in data[colonne_2 - 1] if str(y)!="nan"]
	if colonne_1 == "4":
		XX = [ float('%.2f' % elem) for elem in X ]
	else:
		XX = [ float('%.8f' % elem) for elem in X ]
	for i in range(len(XX)):
		dict_Y[XX[i]] = []
		dict_Y[XX[i]].append(Y[i])
	for i in dict_Y:
		dict_Y[i] = sum(dict_Y[i])/len(dict_Y[i])
	mod_Y = OrderedDict(sorted(dict_Y.items(), key=lambda t: t[0]))

	xy2 = list(zip(mod_Y.keys(), mod_Y.values()))
	poly = xy1 + xy2[::-1]
	lr = LinearRing(poly)
	if lr.is_valid:
		A = Polygon(lr).area
	else:
		A = Polygon(poly).buffer(0).area + Polygon(poly[::-1]).buffer(0).area
	
	f.append(A*10**5)
	print A*10**5

os.chdir("/home/christian/Documents/MyThesis/PyFun")
print "mean :", sum(f)/len(f), "\n", "std  :", numpy.std(f)
