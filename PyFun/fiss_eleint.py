#!/usr/bin/python

import os, numpy, glob, math, collections, linecache, itertools, sys, csv
from Tkinter import *
from tkFileDialog   import askopenfilename 
import tkMessageBox
import matplotlib.pyplot as plt

counting = 0

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def runprog():
    global filename0, filename1, filename2, directory0, directory1, directory2, directory3, filename3, \
           nbr_fiss_sup, nbr_fiss_inf, moy_fiss_sup, moy_fiss_inf, sum_fiss_sup, sum_fiss_inf, esp_fiss_sup, esp_fiss_inf, \
           courbe_nbr_inf, courbe_nbr_sup, courbe_nbr_moy, courbe_ouv_inf, courbe_ouv_sup, courbe_ouv_moy, \
           courbe_mov_inf, courbe_mov_sup, courbe_mov_moy, courbe_esp_inf, courbe_esp_sup, courbe_esp_moy, \
           pas_min, pas_max, pas_int, L_sup, L_inf, seuil, population_fiss

    pas_min = int(ed.get())
    pas_max = int(ef.get())
    pas_int = int(ep.get())
    
    os.chdir(directory0)

    noeuds = int(''.join([k for k in linecache.getline(filename0, 16) if k.isdigit()]))

    os.chdir(directory1)

    coor = {}
    for i in range(noeuds):
        a = ""
        for j in linecache.getline(filename1, i+3):
            if j.isdigit() or j == "." or j == "E" or j == "-":
                a += j
            else:
                if a.isdigit():
                    b = int(a)
                    coor[b]=[]
                elif is_number(a):
                    coor[b].append(float(a))
                a = ""

    surf_sup = max([i[1] for i in coor.values()])
    surf_inf = min([i[1] for i in coor.values()])
                
    noeuds_sup = {}
    noeuds_inf = {}

    for i in [j for j in coor if coor[j][1] == surf_sup]:
        noeuds_sup[i] = coor[i]
    for i in [j for j in coor if coor[j][1] == surf_inf]:
        noeuds_inf[i] = coor[i]

    x_sup_max = max([i[0] for i in noeuds_sup.values()])
    x_sup_min = min([i[0] for i in noeuds_sup.values()])
    x_inf_max = max([i[0] for i in noeuds_inf.values()])
    x_inf_min = min([i[0] for i in noeuds_inf.values()])

    L_sup = x_sup_max - x_sup_min
    L_inf = x_inf_max - x_inf_min

    noeuds_sup_fiss = []
    noeuds_inf_fiss = []

    while len(noeuds_sup) > 0:
        i = noeuds_sup.keys()[0]
        a = []
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

    courbe_nbr_inf = []
    courbe_nbr_sup = []
    courbe_nbr_moy = []
    courbe_ouv_inf = []
    courbe_ouv_sup = []
    courbe_ouv_moy = []
    courbe_mov_inf = []
    courbe_mov_sup = []
    courbe_mov_moy = []
    courbe_esp_inf = []
    courbe_esp_sup = []
    courbe_esp_moy = []

    seuil = float(es.get())

    count = 0
    disp = []
    start = 3 + 3*(pas_max - pas_int) + 4*(pas_max - pas_int) + 2*(pas_max - pas_int)*noeuds              #equation change if it's 3D!!!
    #finish = 3 + 3*(pas_max - pas_int) + 4*(pas_max - pas_int) + (2*(pas_max - pas_int) + 1)*noeuds       #equation change if it's 3D!!!
        
    for i in range(noeuds):
        k = linecache.getline(filename2, i+1+start)
        for j in noeuds_inf_fiss:
            if (i+1) in j:
                disp.append(k)
                break
        for j in noeuds_sup_fiss:
            if (i+1) in j:
                disp.append(k)
                break

    noeuds_disp = {}
    for j in disp:
        b = int(j[:12].replace(" ", ""))
        noeuds_disp[b] = []
        noeuds_disp[b].append(float(j[12:38].replace(" ", "")))
        noeuds_disp[b].append(float(j[38:64].replace(" ", "")))

    for i in noeuds_inf_fiss:
        fiss = 0
        maxi = 0
        j = list(itertools.combinations(i, 2))
        for k in j:
            fiss = math.sqrt((noeuds_disp[k[0]][0] - noeuds_disp[k[1]][0])**2 + (noeuds_disp[k[0]][1] - noeuds_disp[k[1]][1])**2)
            if fiss > maxi:
                maxi = fiss
        if maxi < seuil:
            noeuds_inf_fiss.remove(i)

    for i in noeuds_sup_fiss:
        fiss = 0
        maxi = 0
        j = list(itertools.combinations(i, 2))
        for k in j:
            fiss = math.sqrt((noeuds_disp[k[0]][0] - noeuds_disp[k[1]][0])**2 + (noeuds_disp[k[0]][1] - noeuds_disp[k[1]][1])**2)
            if fiss > maxi:
                maxi = fiss
        if maxi < seuil:
            noeuds_sup_fiss.remove(i)
    
    population_fiss = {}
    
    for z in range(pas_min, pas_max, pas_int):
        disp = []
        start = 3 + 3*z + 4*z + 2*z*noeuds              #equation change if it's 3D!!!
        #finish = 3 + 3*z + 4*z + (2*z + 1)*noeuds       #equation change if it's 3D!!!
        for i in range(noeuds):
            k = linecache.getline(filename2, i+1+start)
            for j in noeuds_inf_fiss:
                if (i+1) in j:
                    disp.append(k)
                    break
            for j in noeuds_sup_fiss:
                if (i+1) in j:
                    disp.append(k)
                    break

        noeuds_disp = {}
        for j in disp:
            b = int(j[:12].replace(" ", ""))
            noeuds_disp[b] = []
            noeuds_disp[b].append(float(j[12:38].replace(" ", "")))
            noeuds_disp[b].append(float(j[38:64].replace(" ", "")))

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
                
        if len(ouv_fiss_sup) > 0:
            moy_fiss_sup = sum(ouv_fiss_sup) / float(len(ouv_fiss_sup))
            sum_fiss_sup = sum(ouv_fiss_sup)
        else:
            moy_fiss_sup = 0
            sum_fiss_sup = 0

        if len(ouv_fiss_inf) > 0:
            moy_fiss_inf = sum(ouv_fiss_inf) / float(len(ouv_fiss_inf))
            sum_fiss_inf = sum(ouv_fiss_inf)
        else:
            moy_fiss_inf = 0
            sum_fiss_inf = 0

        if nbr_fiss_sup > 1:
            esp_fiss_sup = abs((max(coo_fiss_sup) - min(coo_fiss_sup)) / (nbr_fiss_sup - 1))
        else:
            esp_fiss_sup = L_sup

        if nbr_fiss_inf > 1:
            esp_fiss_inf = abs((max(coo_fiss_inf) - min(coo_fiss_inf)) / (nbr_fiss_inf - 1))
        else:
            esp_fiss_inf = L_inf

        courbe_nbr_inf.append(nbr_fiss_inf)
        courbe_nbr_sup.append(nbr_fiss_sup)
        courbe_nbr_moy.append((nbr_fiss_sup + nbr_fiss_inf)*0.5)
        courbe_ouv_inf.append(sum_fiss_inf)
        courbe_ouv_sup.append(sum_fiss_inf)
        courbe_ouv_moy.append((sum_fiss_sup + sum_fiss_inf)*0.5)
        courbe_mov_inf.append(moy_fiss_inf)
        courbe_mov_sup.append(moy_fiss_sup)
        courbe_mov_moy.append((moy_fiss_sup + moy_fiss_inf)*0.5)
        courbe_esp_inf.append(esp_fiss_inf)
        courbe_esp_sup.append(esp_fiss_sup)
        if nbr_fiss_inf > 1 and nbr_fiss_inf <= 1:
            courbe_esp_moy.append(esp_fiss_inf)
        elif nbr_fiss_inf >= 1 and nbr_fiss_inf < 1:
            courbe_esp_moy.append(esp_fiss_sup)
        else:
            courbe_esp_moy.append((esp_fiss_sup + esp_fiss_inf)*0.5)
        
        population_fiss[f3[z, 4]] = ouv_fiss_inf			# 3 is for the Force column in .reac... so 4 if poutre
    linecache.clearcache()
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

def callback3():
    global directory3, filename3, f3
    name = askopenfilename()
    directory3 = name
    filename3 = ""
    for i in range(len(name)):
        j = name[len(name) - i - 1]
        if j == '/':
            directory3 = directory3[:-1]
            break
        else:
            filename3 += j
            directory3 = directory3[:-1]
    filename3 = filename3[::-1] 

    er.delete(0, END)
    er.insert(0, filename3)

    os.chdir(directory3)
    f3 = numpy.loadtxt(filename3)
    lenmax = 0
    for i in f3:
        if str(i[0]) == "nan":
            break
        else:
            lenmax += 1

    ed.delete(0, END)
    ed.insert(0, "0")
    ef.delete(0, END)
    ef.insert(0, str(lenmax))

def sel():
    if str(var1.get())=="1":
        n_fiss = nbr_fiss_sup
        #o_fiss = moy_fiss_sup
        o_fiss = sum_fiss_sup
        e_fiss = esp_fiss_sup
    elif str(var1.get())=="2":
        n_fiss = nbr_fiss_inf
        #o_fiss = moy_fiss_inf
        o_fiss = sum_fiss_inf
        e_fiss = esp_fiss_inf
    elif str(var1.get())=="3":
        n_fiss = (nbr_fiss_sup + nbr_fiss_inf)/2
        #o_fiss = (moy_fiss_sup + moy_fiss_inf)/2
        o_fiss = (sum_fiss_sup + sum_fiss_inf)/2
        e_fiss = (esp_fiss_sup + esp_fiss_inf)/2

    e3.delete(0, END)
    e3.insert(0, str(n_fiss))
    e4.delete(0, END)
    e4.insert(0, str(o_fiss))
    e5.delete(0, END)
    e5.insert(0, str(e_fiss))

def draw():
    global counting
    counting += 1

    counter = 0
    os.chdir(directory3)
    colonne_y = int(ec.get())
    fichier = glob.glob(filename3)[0]
    datalist = numpy.loadtxt(fichier)
    n = len(datalist[0])
    data = []
    for i in range(n):
        data.append(datalist[:,i])
    Y = []
    for z in range(pas_min, pas_max, pas_int):
        Y.append(data[colonne_y - 1][z])
        
    if str(var1.get())=="1":
        X1 = courbe_nbr_sup
        X2 = courbe_ouv_sup
        X3 = courbe_esp_sup
        X4 = courbe_mov_sup
        abscisse1 = "nombre de fissures, face superieur"
        abscisse2 = "ouverture des fissures, face superieur"
        abscisse3 = "espacement des fissures, face superieur"
        abscisse4 = "ouv moy des fissures, face superieur"
    elif str(var1.get())=="2":
        X1 = courbe_nbr_inf
        X2 = courbe_ouv_inf
        X3 = courbe_esp_inf
        X4 = courbe_mov_inf
        abscisse1 = "nombre de fissures, face inferieur"
        abscisse2 = "ouverture des fissures, face inferieur"
        abscisse3 = "espacement des fissures, face inferieur"
        abscisse4 = "ouv moy des fissures, face inferieur"
    elif str(var1.get())=="3":
        X1 = courbe_nbr_moy
        X2 = courbe_ouv_moy
        X3 = courbe_esp_moy
        X4 = courbe_mov_moy
        abscisse1 = "nombre de fissures"
        abscisse2 = "ouverture des fissures"
        abscisse3 = "espacement des fissures"
        abscisse4 = "ouv moy des fissures"

    ordonne = str(eo.get())

    f1 = plt.figure(1)
    ax1 = f1.add_subplot(111)
    ax1.plot(Y, X1, label=filename0[:-5] + " ( >" + str(seuil) + " m)")
    plt.xlabel(ordonne)
    plt.ylabel(abscisse1)
    plt.title(abscisse1 + "/" + ordonne)
    plt.grid(True)
    plt.legend(fancybox=True).get_frame().set_alpha(0.5)
    #counter += 1
    #f1.savefig('/home/cnader/Desktop/point-cet-aprem/' + filename0[:-9] + '-(' + str(counting) + '.' + str(counter) + ').png')
    #plt.clf()
    
    f2 = plt.figure(2)
    ax2 = f2.add_subplot(111)
    ax2.plot(Y, X2, label=filename0[:-5] + " ( >" + str(seuil) + " m)")
    plt.xlabel(ordonne)
    plt.ylabel(abscisse2)
    plt.title(abscisse2 + "/" + ordonne)
    plt.grid(True)
    plt.legend(fancybox=True).get_frame().set_alpha(0.5)
    #counter += 1
    #f2.savefig('/home/cnader/Desktop/point-cet-aprem/' + filename0[:-9] + '-(' + str(counting) + '.' + str(counter) + ').png')
    #plt.clf()
    
    f3 = plt.figure(3)
    ax3 = f3.add_subplot(111)
    ax3.plot(Y, X3, label=filename0[:-5] + " ( >" + str(seuil) + " m)")
    plt.xlabel(ordonne)
    plt.ylabel(abscisse3)
    plt.title(abscisse3 + "/" + ordonne)
    plt.grid(True)
    plt.legend(fancybox=True).get_frame().set_alpha(0.5)
    #counter += 1
    #f3.savefig('/home/cnader/Desktop/point-cet-aprem/' + filename0[:-9] + '-(' + str(counting) + '.' + str(counter) + ').png')
    #plt.clf()
    
    f4 = plt.figure(4)
    ax4 = f4.add_subplot(111)
    ax4.plot(Y, X4, label=filename0[:-5] + " ( >" + str(seuil) + " m)")
    plt.xlabel(ordonne)
    plt.ylabel(abscisse4)
    plt.title(abscisse4 + "/" + ordonne)
    plt.grid(True)
    plt.legend(fancybox=True).get_frame().set_alpha(0.5)
    #counter += 1
    #f3.savefig('/home/cnader/Desktop/point-cet-aprem/' + filename0[:-9] + '-(' + str(counting) + '.' + str(counter) + ').png')
    #plt.clf()
    
    plt.show()
    #plt.close()

def save():
    counter = 0
    for i in directory0[::-1]:
        counter += 1
        if i == "/":
            break

    os.chdir(directory3)
    datalist = numpy.loadtxt(filename3)
    
    os.chdir(directory0[:(len(directory0)-counter)])
    
    if str(var1.get())=="1":
        X1 = courbe_nbr_sup
        X2 = courbe_ouv_sup
        X3 = courbe_esp_sup
        X4 = courbe_mov_sup
    elif str(var1.get())=="2":
        X1 = courbe_nbr_inf
        X2 = courbe_ouv_inf
        X3 = courbe_esp_inf
        X4 = courbe_mov_inf
    else :
        X1 = courbe_nbr_moy
        X2 = courbe_ouv_moy
        X3 = courbe_esp_moy
        X4 = courbe_mov_moy
    
    fn = open(filename0[:-5] + "-" + str(seuil) + "-n-" + directory0[-1] + ".txt", "w")
    counter = pas_min
    for i in X1:
        for j in datalist[counter]:
            fn.write('%.17E' % j + "  ")
        fn.write('%.17E' % i + "\n")
        counter += pas_int
    fn.close()
    fo = open(filename0[:-5] + "-" + str(seuil) + "-o-" + directory0[-1] + ".txt", "w")
    counter = pas_min
    for i in X2:
        for j in datalist[counter]:
            fo.write('%.17E' % j + "  ")
        fo.write('%.17E' % i + "\n")
        counter += pas_int
    fo.close()
    fe = open(filename0[:-5] + "-" + str(seuil) + "-e-" + directory0[-1] + ".txt", "w")
    counter = pas_min
    for i in X3:
        for j in datalist[counter]:
            fe.write('%.17E' % j + "  ")
        fe.write('%.17E' % i + "\n")
        counter += pas_int
    fe.close()
    fom = open(filename0[:-5] + "-" + str(seuil) + "-om-" + directory0[-1] + ".txt", "w")
    counter = pas_min
    for i in X4:
        for j in datalist[counter]:
            fom.write('%.17E' % j + "  ")
        fom.write('%.17E' % i + "\n")
        counter += pas_int
    fom.close()
    fp = csv.writer(open(filename0[:-5] + "-population_fiss.csv", 'wb'))
    for key in population_fiss.keys():
		fp.writerow([key, population_fiss[key]])

###############################################################################
    
root = Tk()
root.title("Depouillement fissuration")

Button(root, text='Fichier .data', command=callback0).grid(row=0, column=0, sticky=W)
Button(root, text='Fichier .msh', command=callback1).grid(row=1, column=0, sticky=W)
Button(root, text='Fichier .res  ', command=callback2).grid(row=2, column=0, sticky=W)
Button(root, text='Fichier .reac', command=callback3).grid(row=3, column=0, sticky=W)

Label(root, text="Seuil").grid(row=4, column=0)
Label(root, text="Pas de calcul:").grid(row=5)
Label(root, text="Nbr de fissures").grid(row=8, column=0)
Label(root, text="Ouv de fissures").grid(row=9, column=0)
Label(root, text="Esp de fissures").grid(row=10, column=0)
Label(root, text="Colonne n:").grid(row=11, column=0, sticky=W)

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
er = Entry(root, width = 26)
es = Entry(root, width = 10, justify = CENTER)
e3 = Entry(root, width = 20, justify = CENTER)
e4 = Entry(root, width = 20, justify = CENTER)
e5 = Entry(root, width = 20, justify = CENTER)
ed = Entry(root, width = 5, justify = CENTER)
ef = Entry(root, width = 5, justify = CENTER)
ep = Entry(root, width = 5, justify = CENTER)
ec = Entry(root, width = 3, justify = CENTER)
eo = Entry(root, width = 20, justify = CENTER)

var1 = IntVar()
R1 = Radiobutton(root, text="sup", variable=var1, value=1, command=sel).grid(row=7, column=0)
R2 = Radiobutton(root, text="inf", variable=var1, value=2, command=sel).grid(row=7, column=1, sticky=W)
R3 = Radiobutton(root, text="moy", variable=var1, value=3, command=sel).grid(row=7, column=1)

e0.grid(row=0, column=1)
e1.grid(row=1, column=1)
e2.grid(row=2, column=1)
er.grid(row=3, column=1)
es.grid(row=4, column=1, sticky=W)
e3.grid(row=8, column=1, sticky=W)
e4.grid(row=9, column=1, sticky=W)
e5.grid(row=10, column=1, sticky=W)
ed.grid(row=5, column=1, sticky=W)
ef.grid(row=5, column=1)
ep.grid(row=5, column=1, sticky=E)
ec.grid(row=11, column=0, sticky=E)
eo.grid(row=11, column=1, sticky=W)

e0.insert(0, "")
e1.insert(0, "")
e2.insert(0, "")
er.insert(0, "")
es.insert(0, "3.0e-5")
e3.insert(0, "")
e4.insert(0, "")
e5.insert(0, "")
ed.insert(0, "Debut")
ef.insert(0, "Fin")
ep.insert(0, "Pas")
ec.insert(0, "1")
eo.insert(0, "Force ou Deplacement?")

Button(root, text='Run code', command=runprog).grid(row=6, column=1, sticky=W, pady=4, columnspan=2)
Button(root, text='Save txt', command=save).grid(row=6, column=1, sticky=E, pady=4, columnspan=2)
Button(root, text='trace courbe', command=draw).grid(row=12, column=1, sticky=W, pady=4, columnspan=2)

mainloop()