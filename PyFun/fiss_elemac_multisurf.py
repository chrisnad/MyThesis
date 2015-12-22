import os, numpy, math, linecache, sys
import operator
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

filename = "poutre_macro_last_CN2"
directory = "/home/christian/Documents/FIDES/FIDES_POST/FIDES_resultats/" + filename + "/3"

os.chdir(directory)

filename0 = filename + ".data"
filename1 = filename + ".gid.msh"
filename2 = filename + ".gid.res"
filename3 = filename + ".reac"

pas_min = 0
pas_max = 530
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

surf = [surf_inf, 0.013815, 0.02763, 0.06963, 0.08, 0.09037]

ouv = []
nbr = []
ouvmoy = []
abs_ouv = []

for I in surf:
	noeuds_surf = {}
	noeuds_surf_disp = {}
	noeuds_surf_pos = {}
	delta = {}
	DP = {}

	for i in [j for j in coor if coor[j][1] == I]:
		noeuds_surf[i] = coor[i]

	x_inf_max = max([i[0] for i in noeuds_surf.values()])
	x_inf_min = min([i[0] for i in noeuds_surf.values()])
	L_inf = x_inf_max - x_inf_min

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
			if int(l[0]) in noeuds_surf.keys():
				if int(l[0]) in noeuds_surf_disp.keys():
					noeuds_surf_disp[int(l[0])].append(l[1:])
				else:
					noeuds_surf_disp[int(l[0])] = [l[1:]]

	for i in noeuds_surf_disp:
		noeuds_surf_pos[i] = []
		for j in noeuds_surf_disp[i]:
			noeuds_surf_pos[i].append([x+y for x,y in zip(j, noeuds_surf[i])])
	
	noeuds_surf_pos = sorted(noeuds_surf_pos.items(), key=lambda k: k[1][0])
	
	for i in range(len(noeuds_surf_pos)-1):
		A = noeuds_surf_pos[i][1]
		B = noeuds_surf_pos[i+1][1]
		Lem = math.sqrt((B[0][0]-A[0][0])**2 + (B[0][1]-A[0][1])**2)
		absci = (B[0][0] + A[0][0])/2.0
		delta[absci] = []
		for j in range(len(A)):
			delta[absci].append(math.sqrt((B[j][0]-A[j][0])**2 + (B[j][1]-A[j][1])**2) - Lem)
	
	c = 0
	for z in range(pas_min, pas_max, pas_int):
		d = -reac[z, 4]							# pour poutres
        #d = reac[z, 3]							# pour tirants
		DP[d] = []
		for i in delta.values():
			DP[d].append(i[c])
		c += 1
	
	DP = sorted(DP.items(), key=operator.itemgetter(0))
	delta = sorted(delta.items(), key=operator.itemgetter(0))
	
	AO = []
	for i in range(len(range(pas_min,pas_max,pas_int))):
		AO.append([i])
	
	c = 0
	for z in range(pas_min,pas_max,pas_int):
		d = -reac[z, 4]							# pour poutres
		#d = reac[z, 3]							# pour tirants
		AO[c].append(d)
		e = []
		for i in delta:
			e.append(i[1][c])
		AO[c].append(e)
		c += 1

	ouv_i = []
	inistep = 0.000009
	for n in DP:
		ouv_i.append([n[0], sum([i for i in n[1] if i>inistep]), I])

	nbr_i = []
	for n in DP:
		nbr_i.append([n[0], len([i for i in n[1] if i>inistep]), I])

	ouvmoy_i = []
	for n in DP:
		ouvmoy_i.append([n[0], avg([i for i in n[1] if i>inistep]), I])

	abs_ouv_i = []
	for i in range(len(range(pas_min,pas_max,pas_int))):
		for n in delta:	
			abs_ouv.append([n[0], n[1][i], i+1, I])
			
	for i in ouv_i:
		ouv.append(i)
	for i in nbr_i:
		nbr.append(i)
	for i in ouvmoy_i:
		ouvmoy.append(i)
	for i in abs_ouv_i:
		abs_ouv.append(i)

linecache.clearcache()
os.chdir(directory[:-2])

f1 = plt.figure(1)
ax1 = f1.add_subplot(111)
ax1.plot([i[0] for i in ouv], [i[1] for i in ouv])
fo = open(filename + "-o-" + directory[-1] + ".txt", "w")
fo.write("Force  Ouverture_de_fissures  Ordonnee\n")
for i in ouv:
	fo.write('%.17E' % i[0] + "  ")
	fo.write('%.17E' % i[1] + "  ")
	fo.write('%.17E' % i[2] + "\n")
fo.close()

f2 = plt.figure(2)
ax2 = f2.add_subplot(111)
ax2.plot([i[0] for i in nbr], [i[1] for i in nbr])
fn = open(filename + "-n-" + directory[-1] + ".txt", "w")
fn.write("Force  Nombre_de_fissures  Ordonnee\n")
for i in nbr:
	fn.write('%.17E' % i[0] + "  ")
	fn.write('%.17E' % i[1] + "  ")
	fn.write('%.17E' % i[2] + "\n")
fn.close()

f3 = plt.figure(3)
ax3 = f3.add_subplot(111)
ax3.plot([i[0] for i in ouvmoy], [i[1] for i in ouvmoy])
fom = open(filename + "-om-" + directory[-1] + ".txt", "w")
fom.write("Force  Ouverture_moyenne_de_fissures  Ordonnee\n")
for i in ouvmoy:
	fom.write('%.17E' % i[0] + "  ")
	fom.write('%.17E' % i[1] + "  ")
	fom.write('%.17E' % i[2] + "\n")
fom.close()

f4 = plt.figure(4)
ax4 = f4.add_subplot(111)
ax4.plot([i[0] for i in abs_ouv[-80:]], [i[1] for i in abs_ouv[-80:]])
fao = open(filename + "-ao-" + directory[-1] + ".txt", "w")
fao.write("Abscisse  Ouverture_de_fissures  Pas_de_calcul Force Ordonne\n")
for i in abs_ouv:
	fao.write('%.17E' % i[0] + "  ")
	fao.write('%.17E' % i[1] + "  ")
	fao.write('%.17E' % i[2] + "  ")
	fao.write('%.17E' % AO[i[2]-1][1] + "  ")
	fao.write('%.17E' % i[3] + "\n")
fao.close()

os.chdir("/home/christian/Documents/MyThesis/PyFun")
plt.show()
