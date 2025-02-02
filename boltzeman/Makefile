F90 = gfortran
FLAGS = -O0 -march=native -g -Wall -ffpe-trap=invalid,overflow -fcheck=all -pedantic
EXE = bgk
OBJ = bgk.o mod_parametres.o
$(EXE) : $(OBJ)
	$(F90) $(FLAGS) -o $(EXE) $^
mod_parametres.o : mod_parametres.f90 
	$(F90) $(FLAGS) -c $<
bgk.o : bgk.f90 mod_parametres.o 
	$(F90) $(FLAGS) -c $<

clean :
	rm -f *.o *.mod $(EXE)