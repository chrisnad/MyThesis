#!/usr/bin/python

import numpy, os
from scipy.optimize import fmin
from shapely.geometry.polygon import  LinearRing, Polygon

global Em, dep_max, xy, alpha, gamma

def area(p):
    return 0.5*abs(sum(x0*y1 - x1*y0 for ((x0, y0), (x1,y1)) in segments(p)))

def segments(p):
    return zip(p, p[1:]+[p[0]])

def obj_fun(X):
    beta = X[0]
    delta = X[1]

    gamma = (2*Ae - Em*alpha*beta - dep_max*delta + delta*beta)/(dep_max - alpha)

    xy0 = [[0.0,0.0]]
    xy0.append([alpha, Em*alpha])
    xy0.append([beta,gamma])
    xy0.append([dep_max, delta])
    
    poly = xy + xy0[::-1]

    lr = LinearRing(poly)
    if lr.is_valid:
        A = Polygon(lr).area
    else:
        A = Polygon(poly).buffer(0).area + Polygon(poly[::-1]).buffer(0).area
        
    return A

XA = []
YA = []
XB = []
YB = []
XC = []
YC = []

for J in range (1, 16):
    
    os.chdir("/home/cnader/Desktop/res_global")
    f = numpy.loadtxt("tirant_macro_3.0-" + str(J) + ".reac")
    x = list(numpy.array([i[0] for i in f]))
    y = list(numpy.array([i[3] for i in f]))
    Em = y[1]/x[1]
    Ef = 3366
    dep_max = x[-1]
    xy = list(zip(x, y))
    alpha = 0
    for i in range(len(xy)):
        if abs(((xy[i][1]/xy[i][0]) - Em)/Em) >  0.25:
            alpha = xy[i-1][0]
            break
    Ae = area(xy + [[x[-1],0]])

    x0 = [0.00082, 0.107]
    
    xopt = fmin(obj_fun, x0)
    beta = xopt[0]
    delta = xopt[1]

    gamma = (2*Ae - Em*alpha*beta - dep_max*delta + delta*beta)/(dep_max - alpha)

    XA.append(alpha)
    YA.append(Em*alpha)
    XB.append(beta)
    YB.append(gamma)    
    XC.append(dep_max)
    YC.append(delta)
    
    x1 = list(numpy.linspace(min(x), max(x), num=1000))
    y1 = []
    for xi in x1:
        if xi <= alpha:
            y1.append(Em*xi)
        elif xi <= beta:
            a = (gamma - Em*alpha)/(beta - alpha)
            b = Em*alpha - a*alpha
            y1.append(a*xi + b)
        elif xi <= delta:
            a = (delta - gamma)/(dep_max - beta)
            b = delta - a*dep_max
            y1.append(a*xi + b)
    x1 = x1[:len(y1)]

    os.chdir("/home/cnader/Desktop/courbes_model_glob")

    f = open("modele_tirant_macro_3.0-" + str(J) + ".txt", "w")
    for j in range(len(x1)):
        f.write('%.17E' % x1[j] + "  ")
        f.write('%.17E' % y1[j] + "\n")
    f.close()

f = open("Variables.txt", "w")
for j in range(len(XA)):
    f.write('%.17E' % XA[j] + " ")
    f.write('%.17E' % YA[j] + "    ")
    f.write('%.17E' % XB[j] + " ")
    f.write('%.17E' % YB[j] + "    ")
    f.write('%.17E' % XC[j] + " ")
    f.write('%.17E' % YC[j] + "\n")
f.close()
