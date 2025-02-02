F90 = gfortran
FLAGS = -O0 -march=native -g -Wall -ffpe-trap=invalid,overflow -fcheck=all -pedantic -fopenmp
LDFLAGS = -fopenmp
EXE = bgk
OBJ = bgk.o mod_parametres.o
DATA_FILES = $(wildcard *.dat)  # Récupère tous les fichiers .dat dans le répertoire

$(EXE) : $(OBJ)
	$(F90) $(FLAGS) -o $(EXE) $^ $(LDFLAGS)

mod_parametres.o : mod_parametres.f90
	$(F90) $(FLAGS) -c $<

bgk.o : bgk.f90 mod_parametres.o
	$(F90) $(FLAGS) -c $<

clean :
	rm -f *.o *.mod $(EXE) $(DATA_FILES)  # Supprime les fichiers .o, .mod, l'exécutable et les fichiers .dat