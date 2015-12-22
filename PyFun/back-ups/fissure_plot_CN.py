#!/usr/bin/python

# TRACE COURBES RESULTATS FICHIERS .REAC

import os, numpy, pylab, glob #sudo apt-get install python-matplotlib
from Tkinter import *
from tkFileDialog   import askopenfilename, asksaveasfilename
import tkMessageBox
import matplotlib.pyplot as plt
from collections import *

yes_mean = -1
one_draw = -1

class Fetch:
    def __init__(self, in_put):
        self.Fetch = in_put

    def __str__(self):
        return self.Fetch

    def get_str(self):
        return str(self.Fetch)

    def get_int(self):
        return int(self.Fetch)

    def get_float(self):
        return float(self.Fetch[1:])

def fetch():
    global filename, directory, colonne_1, colonne_2, abscisse, ordonnee, extention, x_mult, y_mult, title
    directory = Fetch(e0.get()).get_str()
    filename = Fetch(e1.get()).get_str()
    colonne_1 = Fetch(e2.get()).get_int()
    colonne_2 = Fetch(e3.get()).get_int()
    abscisse = Fetch(e4.get()).get_str()
    ordonnee = Fetch(e5.get()).get_str()
    extention = Fetch(e6.get()).get_str()
    x_mult = Fetch(e7.get()).get_float()
    y_mult = Fetch(e8.get()).get_float()
    title = Fetch(e9.get()).get_str()

    if e2.get() == "1":
        e4.delete(0, END)
        e4.insert(0, "Deplacement (m)")

    e9.delete(0, END)
    e9.insert(0, e5.get() + " / " + e4.get())

def draw():
    global f
    os.chdir(directory)
    mean_Y = {}
    f = plt.figure(1)
    ax = f.add_subplot(111)
    
    for fichier in glob.glob(filename + "*" + "*" + extention):
        datalist = numpy.loadtxt(fichier)
        n = len(datalist[0])
        data = []
        for i in range(n):
            data.append(datalist[:,i])
        X = [x*x_mult for x in data[colonne_1 - 1] if str(x)!="nan"]
        Y = [y*y_mult for y in data[colonne_2 - 1] if str(y)!="nan"]
        if yes_mean == 1 or one_draw == 1:
            if e2.get() == "4":
                XX = [ float('%.2f' % elem) for elem in X ]
            else:
                XX = [ float('%.8f' % elem) for elem in X ]
            for i in range(len(XX)):
                if XX[i] in mean_Y:
                    mean_Y[XX[i]].append(Y[i])
                else:
                    mean_Y[XX[i]] = []
                    mean_Y[XX[i]].append(Y[i])
        if one_draw == -1:
            ax.plot(X,Y)

    if mean_Y!={}:
        for i in mean_Y:
            mean_Y[i] = sum(mean_Y[i])/len(mean_Y[i])

        mean_Y = OrderedDict(sorted(mean_Y.items(), key=lambda t: t[0]))    
        ax.plot(mean_Y.keys(), mean_Y.values(), '--', label=filename)
        plt.legend(fancybox=True).get_frame().set_alpha(0.5)
        
    plt.xlabel(abscisse)
    plt.ylabel(ordonnee)
    plt.title(title)
    plt.grid(True)
    plt.show()
    os.chdir("/home/cnader/Documents/Python_functions")

def clear():
    plt.cla()
    pylab.show()

def mean_or_no():
    global yes_mean
    yes_mean *= -1

def only_mean():
    global one_draw
    one_draw *= -1

def callback():
    name = askopenfilename() 

    for i in range(len(name)):
        j = len(name)-i-1
        filename = ''
        if name[j] == '/':
            for k in range(j+1, len(name)):
                filename += list(name).pop(k)
            break

    for i in range(len(filename)):
        j = len(filename)-i-1
        extension = ''
        if filename[j] == '.':
            for k in range(j, len(filename)):
                extension += list(filename).pop(k)
            break

    fichier = ''
    for i in range(len(filename)-len(extension)):
        fichier += filename[i]

    direction = ''
    for i in range(len(name)-len(filename)-1):
        direction += name[i]

    e0.delete(0, END)
    e0.insert(0, direction)
    e1.delete(0, END)
    e1.insert(0, fichier)
    e6.delete(0, END)
    e6.insert(0, extension)

    if list(fichier[::-1])[2] == "o":
        e5.delete(0, END)
        e5.insert(0, "Ouverture de fissure (m)")
    elif list(fichier[::-1])[2] == "e":
        e5.delete(0, END)
        e5.insert(0, "Espacement de fissures (m)")
    elif list(fichier[::-1])[2] == "n":
        e5.delete(0, END)
        e5.insert(0, "Nbr de fissures")

def save():
    name = asksaveasfilename()
    nomfichier = ""
    for i in name[::-1]:
        if i != "/":
            nomfichier += i
        else:
            break
    nomfichier = nomfichier[::-1]
    
    os.chdir(name[:(len(name)-len(nomfichier)-1)])
    plt.savefig(nomfichier)
    
#########################################################################
    
root = Tk()
root.title("Trace courbes CN")
Label(root, text="file location").grid(row=0)
Label(root, text="file name / extention").grid(row=1)
Label(root, text="colonne des X dans .(ext)").grid(row=2)
Label(root, text="colonne des Y dans .(ext)").grid(row=3)
Label(root, text="Title:").grid(row=4, column=1, columnspan=2)

menubar = Menu(root)
filemenu = Menu(menubar, tearoff=0)
#filemenu.add_command(label="New", command=donothing)
filemenu.add_command(label="Open", command=callback)
#filemenu.add_command(label="Save", command=donothing)
#filemenu.add_command(label="Save as...", command=donothing)
#filemenu.add_command(label="Close", command=donothing)
filemenu.add_separator()
filemenu.add_command(label="Exit", command=root.quit)
menubar.add_cascade(label="File", menu=filemenu)
#menubar.add_cascade(label="Edit", menu=editmenu)
helpmenu = Menu(menubar, tearoff=0)
#helpmenu.add_command(label="Help Index", command=donothing)
#helpmenu.add_command(label="About...", command=donothing)
menubar.add_cascade(label="Help", menu=helpmenu)
root.config(menu=menubar)

e0 = Entry(root, width = 26)
e1 = Entry(root, width = 15)
e2 = Entry(root, width = 4, justify = CENTER)
e3 = Entry(root, width = 4, justify = CENTER)
e4 = Entry(root, width = 14)
e5 = Entry(root, width = 14)
e6 = Entry(root, width = 10)
e7 = Entry(root, width = 6, justify = CENTER)
e8 = Entry(root, width = 6, justify = CENTER)
e9 = Entry(root, width = 16, justify = CENTER)

cb1 = Checkbutton(root, text="draw batchs mean curve", command=mean_or_no).grid(row=4, sticky=N)
cb2 = Checkbutton(root, text="draw medium curve only", command=only_mean).grid(row=5, sticky=S)

e0.grid(row=0, column=1, columnspan=2, sticky=W)
e1.grid(row=1, column=1, columnspan=2, sticky=W)
e2.grid(row=2, column=1, sticky=W)
e3.grid(row=3, column=1, sticky=W)
e4.grid(row=2, column=2, sticky=W)
e5.grid(row=3, column=2, sticky=W)
e6.grid(row=1, column=1, columnspan=2, sticky=E)
e7.grid(row=2, column=1, sticky=E)
e8.grid(row=3, column=1, sticky=E)
e9.grid(row=5, column=1, columnspan=2)

e0.insert(0, "/home/cnader/Documents/FIDES/FIDES_Resultats_calculs/Identification_sur_Tirants/tirant_rond_song")
e1.insert(0, "tirant_rond_song")
e2.insert(0, "4")
e3.insert(0, "8")
e4.insert(0, "Force (MN)")
e5.insert(0, "Nbr de fissures")
e6.insert(0, ".txt")
e7.insert(0, "x1")
e8.insert(0, "x1")
e9.insert(0, "Curvy is okay")

Button(root, text='Get Inputs', command=fetch).grid(row=6, column=1, sticky=W, pady=4)
Button(root, text='Draw', command=draw).grid(row=6, column=0, sticky=E, pady=4)
Button(root, text='Clear', command=clear).grid(row=6, column=0, sticky=W, pady=4)
Button(root, text='Save', command=save).grid(row=6, column=2, sticky=W, pady=4)

mainloop()

###
