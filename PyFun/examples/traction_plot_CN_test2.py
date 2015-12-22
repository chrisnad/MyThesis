# TRACE COURBES RESULTATS FICHIERS .REAC

import os, numpy, pylab, glob #sudo apt-get install python-matplotlib
from Tkinter import *
import tkMessageBox

def fetch():
    global filename, directory, colonne_1, colonne_2, abscisse, ordonnee, extention, x_mult, y_mult
    directory = e0.get()
    filename = e1.get()
    colonne_1 = int(e2.get())
    colonne_2 = int(e3.get())
    abscisse = e4.get()
    ordonnee = e5.get()
    extention = e6.get()
    x_mult = int(e7.get()[1:])
    y_mult = int(e8.get()[1:])

def draw():
    os.chdir(directory)
    if len(glob.glob(filename + "*" + "*" + extention)) > 1 or extention == ".reac":
        for j in range(len(glob.glob(filename + "*" + "*" + extention))):
            fichier = ""
            fichier = filename + str(-j-1) + extention
            datalist = numpy.loadtxt(fichier)
            n = len(datalist[0])
            data = []
            for i in range(n):
                data.append(datalist[:,i])
            X = [x*x_mult for x in data[colonne_1 - 1]]
            Y = [y*y_mult for y in data[colonne_2 - 1]]
            pylab.plot(X, Y)
    else :
        fichier = ""
        fichier = filename + extention
        datalist = numpy.loadtxt(fichier)
        n = len(datalist[0])
        data = []
        for i in range(n):
            data.append(datalist[:,i])
        print data
        pylab.plot(data[colonne_1 - 1], data[colonne_2 -1])
        
    pylab.xlabel(abscisse)
    pylab.ylabel(ordonnee)
    pylab.title(filename)
    pylab.show()

root = Tk()
root.title("Trace courbes CN")
Label(root, text="file location").grid(row=0)
Label(root, text="file name / extention)").grid(row=1)
Label(root, text="colonne des X dans .reac").grid(row=2)
Label(root, text="colonne des Y dans .reac").grid(row=3)

e0 = Entry(root, width = 26)
e1 = Entry(root, width = 15)
e2 = Entry(root, width = 4, justify = CENTER)
e3 = Entry(root, width = 4, justify = CENTER)
e4 = Entry(root, width = 14)
e5 = Entry(root, width = 14)
e6 = Entry(root, width = 10)
e7 = Entry(root, width = 6, justify = CENTER)
e8 = Entry(root, width = 6, justify = CENTER)

e0.grid(row=0, column=1, columnspan=2, sticky=W)
e1.grid(row=1, column=1, columnspan=2, sticky=W)
e2.grid(row=2, column=1, sticky=W)
e3.grid(row=3, column=1, sticky=W)
e4.grid(row=2, column=2, sticky=W)
e5.grid(row=3, column=2, sticky=W)
e6.grid(row=1, column=1, columnspan=2, sticky=E)
e7.grid(row=2, column=1, sticky=E)
e8.grid(row=3, column=1, sticky=E)

e0.insert(0, "/home/cnader/Documents/FIDES/FIDES_Resultats_calculs/ElementdeBase/tirant_HA_1.4.2_allong")
e1.insert(0, "tirant_HA12_1.4.2_allong")
e2.insert(0, "1")
e3.insert(0, "4")
e4.insert(0, "Allongement (mm)")
e5.insert(0, "Force (N)")
e6.insert(0, ".reac")
e7.insert(0, "x1")
e8.insert(0, "x1")

Button(root, text='Get inputs', command=fetch).grid(row=4, column=1, sticky=W, pady=4)
Button(root, text='Draw', command=draw).grid(row=4, column=0, sticky=E, pady=4)


mainloop()

###
