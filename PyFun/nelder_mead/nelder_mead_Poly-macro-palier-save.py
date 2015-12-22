#!/usr/bin/python

import numpy, os
from scipy.optimize import fmin
from shapely.geometry.polygon import  LinearRing, Polygon

global Em, sf, xy, alpha, gamma, delta

def area(p):
    return 0.5*abs(sum(x0*y1 - x1*y0 for ((x0, y0), (x1,y1)) in segments(p)))

def segments(p):
    return zip(p, p[1:]+[p[0]])

def obj_fun(X):
    beta = X[0]
    
    gamma = Em*alpha
    #delta = X[1]
    
    delta = (2*Ae - Em*alpha*beta + alpha*gamma + sf*beta)/(sf + gamma)
    #gamma = (2*Ae - Em*alpha*beta - sf*delta + sf*beta)/(delta - alpha)
    
    xy0 = [[0.0,0.0]]
    xy0.append([alpha, Em*alpha])
    xy0.append([beta,gamma])
    xy0.append([delta,sf])
    
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

os.chdir("/home/christian/Desktop/coupes_tirant/sig-eps")

for filename in os.listdir(os.getcwd()):

    os.chdir("/home/christian/Desktop/coupes_tirant/sig-eps")
    
    f = numpy.loadtxt(filename)
    x = list(numpy.array([i[0] for i in f]))
    y = list(numpy.array([i[1] for i in f]))
    Em = y[1]/x[1]
    Ef = 3366
    sf = y[-1]
    xy = list(zip(x, y))
    alpha = 0
    for i in range(len(xy)):
        if abs(((xy[i][1]/xy[i][0]) - Em)/Em) >  0.1:
            alpha = xy[i-1][0]
            break
    Ae = area(xy + [[x[-1],0]])

    #x0 = [0.00081, 0.01]
    x0 = [0.00081]
    
    xopt = fmin(obj_fun, x0)
    beta = xopt[0]
    
    gamma = Em*alpha
    #delta = xopt[1]
    
    delta = (2*Ae - Em*alpha*beta + alpha*gamma + sf*beta)/(sf + gamma)
    #gamma = (2*Ae - Em*alpha*beta - sf*delta + sf*beta)/(delta - alpha)

    if sf/delta > gamma/beta:
        continue
    
    XA.append(alpha)
    YA.append(Em*alpha)
    XB.append(beta)
    YB.append(gamma)
    XC.append(delta)
    YC.append(sf)
    
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
            a = (sf - gamma)/(delta - beta)
            b = sf - a*delta
            y1.append(a*xi + b)
    x1 = x1[:len(y1)]

    os.chdir("/home/christian/Desktop/courbes_model-palier")

    f = open("modele_palier_" + filename, "w")
    for j in range(len(x1)):
        f.write('%.17E' % x1[j] + "  ")
        f.write('%.17E' % y1[j] + "\n")
    f.close()

f = open("Variables-palier.txt", "w")
for j in range(len(XA)):
    f.write('%.17E' % XA[j] + " ")
    f.write('%.17E' % YA[j] + "    ")
    f.write('%.17E' % XB[j] + " ")
    f.write('%.17E' % YB[j] + "    ")
    f.write('%.17E' % XC[j] + " ")
    f.write('%.17E' % YC[j] + "\n")
f.close()
