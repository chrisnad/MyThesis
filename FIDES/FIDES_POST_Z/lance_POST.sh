#!/bin/bash
# ---------------------------------------------------
#   PREPARATION DU FICHIER DE DONNEES
# ---------------------------------------------------
#
#
# ---------------------------------------------------
#   LANCEMENT DE FIDES
# ---------------------------------------------------
#

# ---------------------------------------------------

# Nbr de calcules à lancer

# Nom du fichier ex : nom.data => nomfich = nom
nomfich='cube_tet'
nomdirsauv='cube_tet'
# Prétraitement = 0, Post-traitement = 1
iPost=1
# Nombre de fichiers pour post-traité
nbrfich=1
# Indice fichier de debut
inddeb=0
#----------------------------------------------------

for ((nfich=$inddeb; nfich<$nbrfich+inddeb; ++nfich )) ;
do

	nomdirstock='FIDES_resultats'
	nomtrav='FIDES_fichiers-sources'

	cd $nomtrav

	#make clean
	make

	cd ..
	#./read.sh $nomfich
	sed '2s/.*/'$nomfich".data"'/' lecture_info.data > tmpfile.data  
	mv tmpfile.data lecture_info.data

	sed '12s/.*/'$iPost'/' lecture_info.data > tmpfile.data 
	mv tmpfile.data lecture_info.data

	sed '16s/.*/'$nfich'/' lecture_info.data > tmpfile.data 
	mv tmpfile.data lecture_info.data

	cd $nomtrav

	./FIDES

	mkdir ../$nomdirstock/$nomdirsauv/

	if [ $iPost == 0 ]
	then
		cp ../$nomdirstock/$nomfich.reac ../$nomdirstock/$nomdirsauv/
		cp ../$nomdirstock/$nomfich.resu ../$nomdirstock/$nomdirsauv/
		cp ../$nomdirstock/$nomfich.list ../$nomdirstock/$nomdirsauv/
	else
		mv ../$nomdirstock/$nomfich.gid.* ../$nomdirstock/$nomdirsauv/$nfich/
		mv ../$nomdirstock/$nomfich*.gid.* ../$nomdirstock/$nomdirsauv/$nfich/
		cp ../../FIDES_data/$nomfich.data ../$nomdirstock/$nomdirsauv/$nfich/
	fi

done
