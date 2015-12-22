#!/usr/bin/python

import itertools, numpy, os
from scipy.optimize import fmin
import matplotlib.pyplot as plt
from shapely.geometry.polygon import  LinearRing, Polygon

global Em, sf, xf, xy, poly, xy0, lr, alpha, gamma

def area(p):
    return 0.5*abs(sum(x0*y1 - x1*y0 for ((x0, y0), (x1,y1)) in segments(p)))

def segments(p):
    return zip(p, p[1:]+[p[0]])
    
os.chdir("/home/christian/Desktop")
f = numpy.loadtxt("ele_macro_trilinear-1.reac")
os.chdir("/home/christian/Documents/MyThesis/PyFun")
aire = 0.1*0.1   # sigma/epsilon instead of force/depl
lem = 0.04125	 # sigma/epsilon instead of force/depl
x = list(numpy.array([i[0]/lem for i in f]))
y = list(numpy.array([i[3]/aire for i in f]))
Em = y[1]/x[1]
sf = y[-1]
xf = x[-1]
xy = list(zip(x, y))
alpha = 0
for i in range(len(xy)):
    if abs(((xy[i][1]/xy[i][0]) - Em)/Em) > 0.2:
        alpha = xy[i-1][0]
        break
Ae = area(xy + [[x[-1],0]])

def obj_fun(X):
    ##alpha = X[0]
    ##beta = X[1]
    ##gamma = X[2]
    ##delta = X[3]
    
    beta = X[0]
    delta = X[1]
    
    ##delta = (2*Ae - Em*alpha*beta + alpha*gamma + sf*beta)/(sf + gamma)
    gamma = (2*Ae - Em*alpha*beta - sf*delta + sf*beta)/(delta - alpha)
    ##beta = (2*Ae + alpha*gamma - sf*delta - gamma*delta)/(Em*alpha - sf)
    ##alpha = (2*Ae + beta*sf - gamma*delta - sf*delta)/(Em*beta - gamma)
    
    xy0 = [[0.0,0.0]]
    xy0.append([alpha, Em*alpha])
    xy0.append([beta,gamma])
    xy0.append([delta,sf])
    xy0.append([xf, sf])
    
    poly = xy + xy0[::-1]

    lr = LinearRing(poly)
    if lr.is_valid:
        A = Polygon(lr).area
    else:
        A = Polygon(poly).buffer(0).area + Polygon(poly[::-1]).buffer(0).area
        
    return A

##x0 = [0.0001, 0.0004, 3, 0.0058]
#x0 = [0.0004, 0.0058]
x0 = [0.0001, 0.001]

xopt = fmin(obj_fun, x0)
beta = xopt[0]
delta = xopt[1]

gamma = (2*Ae - Em*alpha*beta - sf*delta+ sf*beta)/(delta - alpha)
##xopt[3] = (2*Ae - Em*xopt[0]*xopt[1] + xopt[0]*xopt[2] + sf*xopt[1])/(sf + xopt[2])
##xopt[2] = (2*Ae - Em*xopt[0]*xopt[1] - sf*xopt[3] + sf*xopt[1])/(xopt[3] - xopt[0])
##xopt[1] = (2*Ae + xopt[0]*xopt[2] - sf*xopt[3] - xopt[2]*xopt[3])/(Em*xopt[0] - sf)
##xopt[0] = (2*Ae - xopt[3]*xopt[2]+sf*xopt[1]-sf*xopt[3])/(Em*xopt[1] - xopt[2])

def f(X):
    global x1, y1
    beta = X[0]
    delta = X[1]
    #alpha = X[0]
    #beta = X[1]
    #gamma = X[2]
    #delta = X[3]
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
    xy1 = list(zip(x1, y1))
    Ac = area(xy1 + [[x1[-1],0]])
    print abs(Ae - Ac)/Ae, x1[-1], y1[-1], delta, sf
    plt.plot(x1, y1, list(x), list(y))
    plt.show()

f(xopt)
