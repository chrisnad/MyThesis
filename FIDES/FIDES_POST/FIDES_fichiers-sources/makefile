#----------------------------------------------------------------------#
#		   	                   FIDES makefile		   #
#----------------------------------------------------------------------#
F90 = gfortran

OPT = -ffree-form -O4 -c
#-g -fbounds-check -Wall -fbacktrace -finit-real=nan

LIB = /home/christian/Documents/FIDES/FIDES_DEVEL/FIDES_lib

#Include
INCLUDE = -I$(LIB)/include/

#Librairie
MUMPS =-L$(LIB) -ldmumps -lmumps_common -lpord -lmpiseq
LAPACK=$(LIB)/liblapack.a
BLAS=$(LIB)/blas_LINUX.a

OBJ = math.o       \
    formatCSRC.o        \
    variables.o	  \
    sparse_matrix.o	   \
    initialisation.o \
    mumps.o \
    lib_elem.o \
	Maillage.o \
	utilitaire.o \
	gnufor2.o \
	proprietes_aleatoires.o \
	element_interface.o \
	fissuration.o \
	fissuration_BA.o \
	aiguillage.o \
	contraintes.o	\
	calc_elem.o \
	assemblage.o \
	AssemblageCSRC.o \
	maj_contraintes.o \
	efforts_int.o \
    PROB_reducmail.o  \
	post_traitement.o   \
	FIDES.o      \

all : FIDES

FIDES : $(OBJ)
	$(F90) $(OBJ) $(MUMPS) $(LAPACK) $(BLAS) -g -lrt -o FIDES

clean :
	rm -f *.out *.o *.mod FIDES

math.o : math.f
	$(F90) $(OPT) math.f

formatCSRC.o : formatCSRC.f
	$(F90) $(OPT) formatCSRC.f

variables.o : variables.f
	$(F90) $(OPT) variables.f

sparse_matrix.o : sparse_matrix.f
	$(F90) $(OPT) sparse_matrix.f

initialisation.o : initialisation.f
	$(F90) $(OPT) initialisation.f

mumps.o : mumps.f
	$(F90) $(OPT) $(INCLUDE) mumps.f

Maillage.o : Maillage.f
	$(F90) $(OPT) Maillage.f

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

fissuration_BA.o : fissuration_BA.f
	$(F90) $(OPT) fissuration_BA.f

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

maj_contraintes.o : maj_contraintes.f
	$(F90) $(OPT) maj_contraintes.f

efforts_int.o : efforts_int.f
	$(F90) $(OPT) efforts_int.f

PROB_reducmail.o : PROB_reducmail.f
	$(F90) $(OPT) PROB_reducmail.f

post_traitement.o : post_traitement.f
	$(F90) $(OPT) post_traitement.f

FIDES.o : FIDES.f
	$(F90) $(OPT) FIDES.f
