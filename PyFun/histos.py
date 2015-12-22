#!/usr/bin/python

import os, numpy, glob, sys
import matplotlib.pyplot as plt
from shapely.geometry.polygon import  LinearRing, Polygon

def area(p):
    return 0.5*abs(sum(x0*y1 - x1*y0 for ((x0, y0), (x1,y1)) in segments(p)))
    
def segments(p):
    return zip(p, p[1:]+[p[0]])

os.chdir("/home/christian/Documents/MyThesis/Projects/Tirant/ScatterModelPoints/courbes_modele")

f = numpy.loadtxt("Variables.txt")

f1 = plt.figure(1)
ax1 = f1.add_subplot(111)
ax1.plot([i[0] for i in f], [i[1] for i in f], 'o',
         [i[2] for i in f], [i[3] for i in f], 'o',
         [i[4] for i in f], [i[5] for i in f], 'o')
plt.title("distribution des points A, B et C")
plt.grid(True)

E = [numpy.trapz([0.0]+[i[1]]+[i[3]]+[i[5]], [0.0]+[i[0]]+[i[2]]+[i[4]]) for i in f]

f2 = plt.figure(2)
ax2 = f2.add_subplot(111)
ax2.hist(E, bins=40)
plt.title("distribution de l'energie")
plt.grid(True)

os.chdir("/home/christian/Documents/MyThesis/Projects/Tirant/ScatterModelPoints/coupes_tirant/sig-eps")

filename = "tirant_macro_3.0_sig-eps" + "*"

Q = []

for fichier in glob.glob(filename):
	datalist = numpy.loadtxt(fichier)
	n = len(datalist)
	data = []
	for i in range(n):
		data.append(datalist[i])
	
	#a = datalist[-1][1]
	#data.append([a/3400.0, a])
	#data.append([a/3400.0, 0.])
	a = datalist[-1][0]
	data.append([a, 0])
	
	lr = LinearRing(data)
	A = Polygon(lr).area
	Q.append(A)
    
f3 = plt.figure(3)
ax3 = f3.add_subplot(111)
ax3.hist(Q, bins=40)
plt.title("distribution reel de l'energie")
plt.grid(True)

f = open("Energie.txt", "w")
f.write('%d' % len(Q) + " ")
f.write('%d' % 1 + "\n")
for j in range(len(Q)):
    f.write('%.17E' % Q[j] + "\n")
f.close()

plt.show()
