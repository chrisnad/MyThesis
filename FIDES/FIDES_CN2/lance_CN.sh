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
nomfich='tirant_HA12_2.0'
# Prétraitement = 0, Post-traitement = 1
iPost=1
# ---------------------------------------------------

nomdirstock='FIDES_resultats'
nomtrav='FIDES_fichiers-sources'
nomdirsauv='Resu_'$nomfich

cd $nomtrav

#make clean
#make

cd ..
#./read.sh $nomfich
sed '2s/.*/'$nomfich".data"'/' lecture_info.data > tmpfile.data  
mv tmpfile.data lecture_info.data

sed '12s/.*/'$iPost'/' lecture_info.data > tmpfile.data 
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
	cp ../$nomdirstock/$nomfich.gid.* ../$nomdirstock/$nomdirsauv/
fi

