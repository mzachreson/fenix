


# compile information
#FORT=mpif90
#FORT=gfortran
#FORT=ifort
#	$ module swap mpi/openmpi_1.6.3_gnu-4.4 mpi/openmpi_1.6.3_intel-13.0.1
FORT=mpif90  
FLAGS=-O3 -ipo -i-static -no-prec-div -vec_report2 -axAVX,SSE4.2
# FLAGS= -check all -traceback -O0 -g
#FLAGS= -O0 -g
#FLAGS=

# some helper variables
CORE=core
BOUNDS=bounds
MAIN=main

# actual targets and psuedotargets
all: fenix10

fenix10: $(CORE).f90 $(BOUNDS).f90 $(MAIN).f90 initialize.f90 move.f90 collide.f90 reloadcells.f90 dataio.f90 inflow.f90

	$(FORT) $(FLAGS) -c $(CORE).f90
	$(FORT) $(FLAGS) -c inflow.f90
	$(FORT) $(FLAGS) -c $(BOUNDS).f90
	$(FORT) $(FLAGS) -c initialize.f90
	$(FORT) $(FLAGS) -c collide.f90
	$(FORT) $(FLAGS) -c move.f90
	$(FORT) $(FLAGS) -c reloadcells.f90
	$(FORT) $(FLAGS) -c dataio.f90
	$(FORT) $(FLAGS) -c $(MAIN).f90
	$(FORT) $(FLAGS) -o fenix10 $(CORE).o $(BOUNDS).o initialize.o move.o collide.o $(MAIN).o reloadcells.o dataio.o inflow.o
	rm -fr *.mod
	rm -fr *.o

clean:
	rm *.mod
	rm *.o
	rm fenix10





