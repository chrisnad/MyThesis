import os, numpy, pylab, glob, math, matplotlib #sudo apt-get install python-matplotlib
from Tkinter import *
from tkFileDialog   import askopenfilename 
import tkMessageBox
import itertools

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def runprog():
    global filename0, filename1, filename2, directory0, directory1, directory2, nbr_fiss_sup, nbr_fiss_inf, moy_fiss_sup, moy_fiss_inf, esp_fiss_sup, esp_fiss_inf

    os.chdir(directory0)
    f = open(filename0, 'r')

    count = 0
    for i in f:
        count += 1
        if count == 16:
            a = ''
            for j in i:
                if j.isdigit():
                    a += j
                else:
                    if a.isdigit():
                        noeuds = int(a)
                        break
            break

    os.chdir(directory1)
    f = open(filename1, 'r')

    count = 0
    coor = {}
    for i in f:
        count += 1
        if count > 2 and count < 3 + noeuds:
            a = ''
            for j in i:
                if j.isdigit() or j == "." or j == "E" or j == "-":
                    a += j
                else:
                    if a.isdigit():
                        b = int(a)
                        coor[b]=[]
                    elif is_number(a):
                        coor[b].append(float(a))
                    a = ''

    surf_sup = -float("inf")
    surf_inf = float("inf")

    for i in coor:
        if coor[i][1] > surf_sup:
            surf_sup = coor[i][1]
        if coor[i][1] < surf_inf:
            surf_inf = coor[i][1]
                
    noeuds_sup = {}
    noeuds_inf = {}

    for i in coor:
        if coor[i][1] == surf_sup:
            noeuds_sup[i] = coor[i]
        elif coor[i][1] == surf_inf:
            noeuds_inf[i] = coor[i]

    noeuds_sup_fiss = []
    noeuds_inf_fiss = []

    while len(noeuds_sup) > 0:
        i = noeuds_sup.keys()[0]
        a = []
        b = []
        for k in noeuds_sup:
            if k != i:
                if noeuds_sup[i] == noeuds_sup[k]:
                    a.append(i)
                    a.append(k)
        a = list(set(a))
        for j in a:
            del noeuds_sup[j]
        if len(a) != 0:
            noeuds_sup_fiss.append(a)
        else:
            del noeuds_sup[i]

    while len(noeuds_inf) > 0:
        i = noeuds_inf.keys()[0]
        a = []
        for k in noeuds_inf:
            if k != i:
                if noeuds_inf[i] == noeuds_inf[k]:
                    a.append(i)
                    a.append(k)
        a = list(set(a))
        for j in a:
            del noeuds_inf[j]
        if len(a) != 0:
            noeuds_inf_fiss.append(a)
        else:
            del noeuds_inf[i]

    os.chdir(directory2)
    f = open(filename2, 'r')

    count = 0
    disp = []
    for i in f:
        count += 1
        if count > 3:
            for j in noeuds_inf_fiss:
                if (count - 3) in j:
                    disp.append(i)
                    break
            for j in noeuds_sup_fiss:
                if (count - 3) in j:
                    disp.append(i)
                    break

    noeuds_disp = {}
    for i in disp:
        a = ''
        for j in i:
            if j.isdigit() or j == "." or j == "E" or j == "-":
                a += j
            else:
                if a.isdigit():
                    b = int(a)
                    noeuds_disp[b]=[]
                elif is_number(a):
                    noeuds_disp[b].append(float(a))
                a = ''

    ouv_fiss_sup = []
    ouv_fiss_inf = []
    coo_fiss_sup = []
    coo_fiss_inf = []
    nbr_fiss_sup = 0
    nbr_fiss_inf = 0
    seuil = float(es.get())

    for i in noeuds_inf_fiss:
        fiss = 0
        maxi = 0
        j = list(itertools.combinations(i, 2))
        for k in j:
            fiss = math.sqrt((noeuds_disp[k[0]][0] - noeuds_disp[k[1]][0])**2 + (noeuds_disp[k[0]][1] - noeuds_disp[k[1]][1])**2)
            if fiss > maxi:
                maxi = fiss
        if maxi >= seuil:
            ouv_fiss_inf.append(maxi)
            nbr_fiss_inf += 1
            coo_fiss_inf.append(coor[i[0]][0])

    for i in noeuds_sup_fiss:
        fiss = 0
        maxi = 0
        j = list(itertools.combinations(i, 2))
        for k in j:
            fiss = math.sqrt((noeuds_disp[k[0]][0] - noeuds_disp[k[1]][0])**2 + (noeuds_disp[k[0]][1] - noeuds_disp[k[1]][1])**2)
            if fiss > maxi:
                maxi = fiss
        if maxi >= seuil:
            ouv_fiss_sup.append(maxi)
            nbr_fiss_sup += 1
            coo_fiss_sup.append(coor[i[0]][0])

    moy_fiss_sup = sum(ouv_fiss_sup) / float(len(ouv_fiss_sup))
    moy_fiss_inf = sum(ouv_fiss_inf) / float(len(ouv_fiss_inf))
    esp_fiss_sup = abs((max(coo_fiss_sup) - min(coo_fiss_sup)) / (nbr_fiss_sup - 1))
    esp_fiss_inf = abs((max(coo_fiss_inf) - min(coo_fiss_inf)) / (nbr_fiss_inf - 1))

###########################################################################

def callback0():
    global directory0, filename0
    name = askopenfilename()
    directory0 = name
    filename0 = ""
    for i in range(len(name)):
        j = name[len(name) - i - 1]
        if j == '/':
            directory0 = directory0[:-1]
            break
        else:
            filename0 += j
            directory0 = directory0[:-1]
    filename0 = filename0[::-1]

    e0.delete(0, END)
    e0.insert(0, filename0)

def callback1():
    global directory1, filename1
    name = askopenfilename()
    directory1 = name
    filename1 = ""
    for i in range(len(name)):
        j = name[len(name) - i - 1]
        if j == '/':
            directory1 = directory1[:-1]
            break
        else:
            filename1 += j
            directory1 = directory1[:-1]
    filename1 = filename1[::-1]

    e1.delete(0, END)
    e1.insert(0, filename1)

def callback2():
    global directory2, filename2
    name = askopenfilename()
    directory2 = name
    filename2 = ""
    for i in range(len(name)):
        j = name[len(name) - i - 1]
        if j == '/':
            directory2 = directory2[:-1]
            break
        else:
            filename2 += j
            directory2 = directory2[:-1]
    filename2 = filename2[::-1] 

    e2.delete(0, END)
    e2.insert(0, filename2)

def sel():
    
    if str(var.get())=="1":
        n_fiss = nbr_fiss_sup
        o_fiss = moy_fiss_sup
        e_fiss = esp_fiss_sup
    elif str(var.get())=="2":
        n_fiss = nbr_fiss_inf
        o_fiss = moy_fiss_inf
        e_fiss = esp_fiss_inf
    else:
        n_fiss = (nbr_fiss_sup + nbr_fiss_inf)/2
        o_fiss = (moy_fiss_sup + moy_fiss_inf)/2
        e_fiss = (esp_fiss_sup + esp_fiss_inf)/2

    e3.delete(0, END)
    e3.insert(0, str(n_fiss))
    e4.delete(0, END)
    e4.insert(0, str(o_fiss))
    e5.delete(0, END)
    e5.insert(0, str(e_fiss))

###############################################################################
    
root = Tk()
root.title("Depouillement fissuration")

Button(root, text='Fichier .data', command=callback0).grid(row=0, column=0)
Button(root, text='Fichier .msh ', command=callback1).grid(row=1, column=0)
Button(root, text='Fichier .res ', command=callback2).grid(row=2, column=0)

Label(root, text="seuil").grid(row=3, column=0)
Label(root, text="Nbr de fissures").grid(row=6, column=0)
Label(root, text="Ouv de fissures").grid(row=7, column=0)
Label(root, text="Esp de fissures").grid(row=8, column=0)

menubar = Menu(root)
filemenu = Menu(menubar, tearoff=0)
filemenu.add_separator()
filemenu.add_command(label="Exit", command=root.quit)
menubar.add_cascade(label="File", menu=filemenu)

helpmenu = Menu(menubar, tearoff=0)

menubar.add_cascade(label="Help", menu=helpmenu)
root.config(menu=menubar)

e0 = Entry(root, width = 26)
e1 = Entry(root, width = 26)
e2 = Entry(root, width = 26)
es = Entry(root, width = 10, justify = CENTER)
e3 = Entry(root, width = 20, justify = CENTER)
e4 = Entry(root, width = 20, justify = CENTER)
e5 = Entry(root, width = 20, justify = CENTER)

var = IntVar()
R1 = Radiobutton(root, text="sup", variable=var, value=1, command=sel).grid(row=5, column=0)
R2 = Radiobutton(root, text="inf", variable=var, value=2, command=sel).grid(row=5, column=1, sticky=W)
R3 = Radiobutton(root, text="moy", variable=var, value=3, command=sel).grid(row=5, column=1)

e0.grid(row=0, column=1)
e1.grid(row=1, column=1)
e2.grid(row=2, column=1)
es.grid(row=3, column=1, sticky=W)
e3.grid(row=6, column=1, sticky=W)
e4.grid(row=7, column=1, sticky=W)
e5.grid(row=8, column=1, sticky=W)

e0.insert(0, "")
e1.insert(0, "")
e2.insert(0, "")
es.insert(0, "3.0e-5")
e3.insert(0, "")
e4.insert(0, "")
e5.insert(0, "")

Button(root, text='Run code', command=runprog).grid(row=4, column=1, sticky=W, pady=4, columnspan=2)

mainloop()

