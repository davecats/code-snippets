F90 = mpifort
INC = /usr/include/

OBJ = typedef.o mpi_transpose.o transpose_datatypes.o 
OBJNC = typedef.o mpi_transpose.o transpose_packunpack.o
flags =  -fall-intrinsics -ffree-line-length-none  -I$(INC) -cpp  -g -fcheck=all -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow -lfftw3
#flags =  -fall-intrinsics -ffree-line-length-none -std=f2008 -I$(INC) -cpp  -O3

transpose_datatypes: $(OBJ)
	$(F90) $(flags) -o $@ $(OBJ)
transpose_packunpack: $(OBJNC)
	$(F90) $(flags) -o $@ $(OBJNC)
%.o : %.f90
	$(F90) $(flags) -c $<
clean: 
	rm *.mod *.o
