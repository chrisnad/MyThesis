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
# CALCULS PROBABILISTES A LA CHAINE ...
#####################################################
compt=1
nbMaxCalc=1
FILENAME='lecture_info.data'
iPost=0
#####################################################
nomdirstock='FIDES_resultats'
nomtrav='FIDES_fichiers-sources'
fichierfides='FIDES_CN'

reset
more=1
j=1
while (($more == 1))
do
	icompt=$compt
	inbMaxCalc=$nbMaxCalc
	if  (($j == 1))
	then
                nomfich='poutre_macro_CN'
                nomdirsauv='poutre_macro_CN'
        	numdirsup='3D'
		more=0
	fi 

	cd $nomtrav

	#
	#
	# Creation du dossier cible CALTOUT
	# --------------------------------
	#mkdir ../$nomdirstock/$nomdirsauv
	mkdir ../../'FIDES_Resultats_calculs'/$numdirsup/
    mkdir ../../'FIDES_Resultats_calculs'/$numdirsup/$nomdirsauv/
    mkdir ../../'FIDES_POST'/'FIDES_resultats'/$nomdirsauv/
	#
	# Initialisations
	# ---------------
	echo '-------------------------------------------------'
	echo ' Lancement des calculs et stockage des resultats'
	echo '-------------------------------------------------'
	echo ' '
	nbtest=1
	
	#
	while (($nbtest))
	do
	   #
	   # Lancement du calcul FIDES
	   #--------------------------
	   echo " "
	   echo "Lancement du calcul FIDES"
	   echo " "
	   echo "Calcul numero = " $icompt
	   echo " "
	   echo " "
	   #ouvre le fichier lecture_info.data et remplace la 2e ligne pour le nom choisi
	   cd 
	   cd 'Documents'/'FIDES'/$fichierfides/
	   #./read.sh $nomfich
	   sed '2s/.*/'$nomfich".data"'/' lecture_info.data > tmpfile.data  
	   mv tmpfile.data lecture_info.data
	   
	   sed '12s/.*/'$iPost'/' lecture_info.data > tmpfile.data 
	   mv tmpfile.data lecture_info.data
	   cd $nomtrav
	   
	   if (($j == 1))
	   then			
		#make clean
		make
	   fi
	   
	   ./FIDES
	   
	   #
	   # Fin calcul FIDES
	   # -------------------------
	   #
	   #
	   # Menage entre deux calculs consecutifs
	   #--------------------------------------
	   ################################################### ####################################################
	   #if (($icompt == 1))
	   #then
		##cp ../FIDES_data/$nomfich.data ../../'FIDES_POST'/'FIDES_resultats_profils'/$nomdirsauv/$nomfich.data
		#cp ../$nomdirstock/$nomfich.reac ../../'FIDES_POST'/'FIDES_resultats'/$nomdirsauv/$nomfich.reac
		#cp ../$nomdirstock/$nomfich.list ../../'FIDES_POST'/'FIDES_resultats'/$nomdirsauv/$nomfich.list
	   #fi
	   ################################################### ####################################################
	   cp ../../'FIDES_data'/$nomfich.data ../'FIDES_Resultats_calculs'/$numdirsup/$nomdirsauv/$nomfich.data
	   cp ../$nomdirstock/$nomfich.reac ../../'FIDES_Resultats_calculs'/$numdirsup/$nomdirsauv/$nomfich"-"$icompt.reac
	   #rm ../$nomdirstock/$nomfich.list
	   #rm ../$nomdirstock/$nomfich.resu
	   #mv ../$nomdirstock/$nomfich.resu ../$nomdirstock/$nomdirsauv/$nomfich"-"$icompt.resu
	   #cp ../$nomdirstock/$nomfich.list ../$nomdirstock/$nomdirsauv/$nomfich"-"$icompt.list
	   #copy to dossier de résultats centralisés
	   mkdir ../../'FIDES_POST'/'FIDES_resultats'/$nomdirsauv/$icompt/
	   cp ../$nomdirstock/$nomfich.reac ../../'FIDES_POST'/'FIDES_resultats'/$nomdirsauv/$icompt/$nomfich.reac
	   mv ../$nomdirstock/$nomfich.list ../../'FIDES_POST'/'FIDES_resultats'/$nomdirsauv/$icompt/$nomfich.list
	   #cp ../$nomdirstock/$nomfich.list ../../'FIDES_Resultats_calculs'/$nomdirsauv/$nomfich"-"$icompt.list	
	   #mv ../$nomdirstock/$nomfich.reac ../$nomdirstock/$nomdirsauv/$nomfich"-"$icompt.reac
	   
	   
	   # Test sur le nombre de calculs valides
	   #--------------------------------------
	   if  (($icompt >= $compt + $inbMaxCalc -1))
	   then
		  nbtest=0
	   else
		  ((icompt=$icompt+1))
	   fi
	   echo " "
	   echo " "
	   echo 'Fin des calculs de ' $nomfich
	done
	
	# sync automatique du dossier résultats 
	# sur un ordinateur branché au réseau une fois les calculs terminés
	#ping -c 1 -t 1 137.121.81.54
	#if [ $? -eq 0 ]; 
	#then
	#	rsync -r ssh /home/dom_dum/FIDES_Resultats_calculs/ 137.121.81.54:/Volumes/Doctorat/Projet/Modelisation/IFSTTAR/FIDES_Resultats_calculs/
	#fi
	
	j=$j+1
done
	cd ..
	   echo " "
	   echo " "

	echo 'Fin des calculs !'
	   echo " "
	   echo " "

