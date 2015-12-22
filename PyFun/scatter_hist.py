#!/usr/bin/python

import numpy as np
import os, math
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

os.chdir("/home/cnader/Desktop/the_scattering/courbes_modele")

f = np.loadtxt("Variables.txt")

### data
x = [i[0] for i in f]
y = [i[1] for i in f]

nullfmt   = NullFormatter()         ### no labels

### definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left+width+0.02

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]

plt.figure(1)

axScatter = plt.axes(rect_scatter)
plt.grid(True)
plt.xlabel("deformation")
plt.ylabel("contrainte")
axHistx = plt.axes(rect_histx)
axHisty = plt.axes(rect_histy)

### no labels
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

### the scatter plot:
axScatter.scatter(x, y)

### now determine nice limits by hand:
xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )

axScatter.set_xlim( (0.00006, 0.00012) )    # point A
#axScatter.set_xlim( (0.0005, 0.0015) )      # point B
#axScatter.set_xlim( (0.004, 0.008) )        # point C

axHistx.hist(x, bins=300, alpha = 0.2)
axHisty.hist(y, bins=300, orientation='horizontal', alpha = 0.2)

axHistx.set_xlim( axScatter.get_xlim() )
axHisty.set_ylim( axScatter.get_ylim() )

plt.title("dispersion sur A")
plt.show()
