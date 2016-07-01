#----------------------------------------------------------------------#
#		   	     FIDES makefile	                       #
#----------------------------------------------------------------------#

F90 = gfortran

###FOR DEBUGGING PURPOSES
#F90 = gfortran -fcheck=all

OPT = -ffree-form -O4 -c
#-g -fbounds-check -Wall -fbacktrace -finit-real=nan

LIB = /home/christian/Documents/FIDES/FIDES_DEVEL/FIDES_lib

#Include
INCLUDE = -I$(LIB)/include/

#Librairie
MUMPS =-L$(LIB) -ldmumps -lmumps_common -lpord -lmpiseq
LAPACK=$(LIB)/liblapack.a
BLAS=$(LIB)/blas_LINUX.a
	
OBJ =   math.o	\
        formatCSRC.o       \
	variables.o	\
	sparse_matrix.o	\
	initialisation.o	\
	mumps.o	\
	lib_elem.o	\
        Maillage.o             \
	utilitaire.o	\
	gnufor2.o	\
	proprietes_aleatoires.o	\
	element_interface.o     \
	fissuration.o           \
	interface_macro.o     \
        acier.o           \
	aiguillage.o	\
	contraintes.o	\
	calc_elem.o	\
	assemblage.o	\
        AssemblageCSRC.o        \
        PROB_reducmail.o        \
	post_traitement.o        \
	depouil_lance.o          \
	FIDES.o                 \

all : FIDES

FIDES : $(OBJ)
	$(F90) $(OBJ) $(MUMPS) $(LAPACK) $(BLAS) -g -lrt -o FIDES

clean : 
	rm -f *.out *.o *.mod FIDES

math.o : math.f
	$(F90) $(OPT) math.f

variables.o : variables.f
	$(F90) $(OPT) variables.f

formatCSRC.o : formatCSRC.f
	$(F90) $(OPT) formatCSRC.f

sparse_matrix.o : sparse_matrix.f
	$(F90) $(OPT) sparse_matrix.f

Maillage.o : Maillage.f
	$(F90) $(OPT) Maillage.f

initialisation.o : initialisation.f
	$(F90) $(OPT) initialisation.f

mumps.o : mumps.f
	$(F90) $(OPT) $(INCLUDE) mumps.f

lib_elem.o : lib_elem.f
	$(F90) $(OPT) lib_elem.f

utilitaire.o : utilitaire.f
	$(F90) $(OPT) utilitaire.f

gnufor2.o : gnufor2.f
	$(F90) $(OPT) gnufor2.f

proprietes_aleatoires.o : proprietes_aleatoires.f
	$(F90) $(OPT) proprietes_aleatoires.f

element_interface.o : element_interface.f
	$(F90) $(OPT) element_interface.f

fissuration.o : fissuration.f
	$(F90) $(OPT) fissuration.f

interface_macro.o : interface_macro.f
	$(F90) $(OPT) interface_macro.f

acier.o : acier.f
	$(F90) $(OPT) acier.f
	
aiguillage.o : aiguillage.f
	$(F90) $(OPT) aiguillage.f

contraintes.o : contraintes.f
	$(F90) $(OPT) contraintes.f

calc_elem.o : calc_elem.f
	$(F90) $(OPT) calc_elem.f

assemblage.o : assemblage.f
	$(F90) $(OPT) assemblage.f

AssemblageCSRC.o : AssemblageCSRC.f
	$(F90) $(OPT) AssemblageCSRC.f

PROB_reducmail.o : PROB_reducmail.f
	$(F90) $(OPT) PROB_reducmail.f
	
post_traitement.o : post_traitement.f
	$(F90) $(OPT) post_traitement.f

depouil_lance.o : depouil_lance.f
	$(F90) $(OPT) depouil_lance.f

FIDES.o : FIDES.f
	$(F90) $(OPT) FIDES.f
