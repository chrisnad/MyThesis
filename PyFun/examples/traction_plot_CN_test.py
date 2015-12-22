# TRACE COURBES RESULTATS FICHIERS .REAC

import os, numpy, pylab, glob #sudo apt-get install python-matplotlib
from Tkinter import *
import tkMessageBox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

def fetch():
    global filename, directory, colonne_1, colonne_2, abscisse, ordonnee
    directory = e0.get()
    filename = e1.get()
    colonne_1 = int(e2.get())
    colonne_2 = int(e3.get())
    abscisse = e4.get()
    ordonnee = e5.get()

def draw():
    global f
    os.chdir(directory)
    for j in range(len(glob.glob("*.reac"))):
        fichier = ""
        fichier = filename + str(-j-1) +".reac"
        datalist = numpy.loadtxt(fichier)
        n = len(datalist[0])
        data = []
        for i in range(n):
            data.append(datalist[:,i])
        f = pylab.plot(data[colonne_1 - 1], data[colonne_2 -1])

    #f.xlabel(abscisse)
    #f.ylabel(ordonnee)
    #f.title(filename)
    dataplot = FigureCanvasTkAgg(f, master=root)
    dataplot.show()

root = Tk()
root.title("Trace courbes CN")
Label(root, text="file location").grid(row=0)
Label(root, text="file name (no extention)").grid(row=1)
Label(root, text="colonne des X dans .reac").grid(row=2)
Label(root, text="colonne des Y dans .reac").grid(row=3)

e0 = Entry(root, width = 26)
e1 = Entry(root, width = 26)
e2 = Entry(root, width = 4, justify = CENTER)
e3 = Entry(root, width = 4, justify = CENTER)
e4 = Entry(root, width = 14)
e5 = Entry(root, width = 14)

e0.grid(row=0, column=1, columnspan=2, sticky=W)
e1.grid(row=1, column=1, columnspan=2, sticky=W)
e2.grid(row=2, column=1, sticky=W)
e3.grid(row=3, column=1, sticky=W)
e4.grid(row=2, column=1)
e5.grid(row=3, column=1)


e0.insert(0, "/home/cnader/Documents/FIDES/FIDES_Resultats_calculs/ElementdeBase/tirant_HA_1.4.2_allong")
e1.insert(0, "tirant_HA12_1.4.2_allong")
e2.insert(0, "1")
e3.insert(0, "4")
e4.insert(0, "Allongement (mm)")
e5.insert(0, "Force (N)")

Button(root, text='Get inputs', command=fetch).grid(row=4, column=0, sticky=E, pady=4)
Button(root, text='Draw', command=draw).grid(row=4, column=0, sticky=W, pady=4)


c = Canvas(root)
c.pack(side=TOP, fill=BOTH, expand=1)

toolbar = NavigationToolbar2TkAgg(c, root)
toolbar.grid(row=8,column=1, sticky=W)

mainloop()

###
