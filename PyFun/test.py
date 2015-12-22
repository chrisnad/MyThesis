import os, numpy, math, linecache
import operator
from itertools import izip
import matplotlib.pyplot as plt

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
        
def avg(x):
	if len(x) == 0:
		return 0.0
	else:
		return sum(x)/float(len(x))

directory = "/home/christian/Documents/FIDES/FIDES_POST_test/FIDES_resultats/poutre_macro_CN/1"
os.chdir(directory)

filename = "poutre_macro_CN"
filename0 = filename + ".data"
filename1 = filename + ".gid.msh"
filename2 = filename + ".gid.res"
filename3 = filename + ".reac"

pas_min = 0
pas_max = 425
pas_int = 1

reac = numpy.loadtxt(filename3)

noeuds = int(''.join([k for k in linecache.getline(filename0, 16) if k.isdigit()]))

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
    
noeuds_inf_disp = {}
noeuds_sup_disp = {}

step = noeuds*2 + 4 + 3              #equation change if it's 3D!!!

for q in range(pas_min, pas_max, pas_int):
	for r in range(noeuds):
		k = linecache.getline(filename2, 4 + q*step + r)
		l = []
		for elem in k.split():
			try:
				l.append(float(elem))
			except ValueError:
				pass
		if int(l[0]) in noeuds_inf.keys():
			if int(l[0]) in noeuds_inf_disp.keys():
				noeuds_inf_disp[int(l[0])].append(l[1:])
			else:
				noeuds_inf_disp[int(l[0])] = [l[1:]]
		if int(l[0]) in noeuds_sup.keys():
			if int(l[0]) in noeuds_sup_disp.keys():
				noeuds_sup_disp[int(l[0])].append(l[1:])
			else:
				noeuds_sup_disp[int(l[0])] = [l[1:]]

noeuds_inf = sorted(noeuds_inf.items(), key=operator.itemgetter(0))
noeuds_sup = sorted(noeuds_sup.items(), key=operator.itemgetter(0))
noeuds_inf_disp = sorted(noeuds_inf_disp.items(), key=operator.itemgetter(0))
noeuds_sup_disp = sorted(noeuds_sup_disp.items(), key=operator.itemgetter(0))

noeuds_inf_pos = noeuds_inf_disp
noeuds_sup_pos = noeuds_sup_disp

d = 0
for i in noeuds_inf_disp:
	c = 0
	for j in i[1]:
		noeuds_inf_pos[d][1][c] = map(sum, zip(j , noeuds_inf[d][1][:2]))
		c += 1
	d += 1

delta = {}
c = 0
for i in range(len(noeuds_inf_pos)-1):
	A = noeuds_inf_pos[i][1]
	B = noeuds_inf_pos[i+1][1]
	Lem = math.sqrt((B[0][0]-A[0][0])**2 + (B[0][1]-A[0][1])**2)
	delta[c] = []
	for j in range(len(A)):
		delta[c].append(math.sqrt((B[j][0]-A[j][0])**2 + (B[j][1]-A[j][1])**2) - Lem)
	c += 1

DP = {}
c = 0
for z in range(pas_min, pas_max, pas_int):
	d = -reac[z, 4]							# 4 for poutre (and negative), 3 for tirant
	DP[d] = []
	for i in delta.values():
		DP[d].append(i[c])
	c += 1
	
DP = sorted(DP.items(), key=operator.itemgetter(0))

ouv = {}
for n in DP:
	ouv[n[0]] = sum([i for i in n[1] if i>0.000009])
ouv = sorted(ouv.items(), key=operator.itemgetter(0))

nbr = {}
for n in DP:
	nbr[n[0]] = len([i for i in n[1] if i>0.000009])
nbr = sorted(nbr.items(), key=operator.itemgetter(0))

ouvmoy = {}
for n in DP:
	ouvmoy[n[0]] = avg([i for i in n[1] if i>0.000009])
ouvmoy = sorted(ouvmoy.items(), key=operator.itemgetter(0))

os.chdir(directory[:-2])

f1 = plt.figure(1)
ax1 = f1.add_subplot(111)
ax1.plot([i[0] for i in ouv], [i[1] for i in ouv])
fo = open(filename + "-" + "-o-" + directory[-1] + ".txt", "w")
fo.write("Force  Ouverture_de_fissures\n")
for i in ouv:
	fo.write('%.17E' % i[0] + "  ")
	fo.write('%.17E' % i[1] + "\n")
fo.close()

f2 = plt.figure(2)
ax2 = f2.add_subplot(111)
ax2.plot([i[0] for i in nbr], [i[1] for i in nbr])
fn = open(filename + "-" + "-n-" + directory[-1] + ".txt", "w")
fn.write("Force  Nombre_de_fissures\n")
for i in nbr:
	fn.write('%.17E' % i[0] + "  ")
	fn.write('%.17E' % i[1] + "\n")
fn.close()

f3 = plt.figure(3)
ax3 = f3.add_subplot(111)
ax3.plot([i[0] for i in ouvmoy], [i[1] for i in ouvmoy])
fom = open(filename + "-" + "-om-" + directory[-1] + ".txt", "w")
fom.write("Force  Ouverture_de_fissures\n")
for i in ouvmoy:
	fom.write('%.17E' % i[0] + "  ")
	fom.write('%.17E' % i[1] + "\n")
fom.close()

plt.show()

os.chdir("/home/christian/Documents/MyThesis/PyFun")
