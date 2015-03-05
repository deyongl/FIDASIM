#System:   64 bit
#Compiler: gfortran or ifort

ifeq ($(FIDASIM_COMPILER),gfortran)
	LFLAGS = -lnetcdff -lnetcdf -lm
	CFLAGS = -Ofast -fopenmp
	DEBUG_CFLAGS = -Og -fbacktrace -fbounds_check -Wall -ffpe-trap=invalid,zero,overflow
endif

ifeq ($(FIDASIM_COMPILER),ifort)
	LFLAGS = -lnetcdff -lnetcdf -limf -lm
	CFLAGS = -O2 -openmp
	DEBUG_CFLAGS = -O0 -g -traceback -debug all -check all -check bounds -fpe:0 -fpstkchk -warn
endif

all: fidasim

debug: CFLAGS = DEBUG_CFLAGS
debug: fidasim_debug

fidasim fidasim_debug: fidasim.o fidasim_linalg.o
	$(FIDASIM_COMPILER) $(CFLAGS) $^ -o $@ -L$(NETCDF_LIB) $(LFLAGS)

fidasim.o: fidasim.f90 fidasim_linalg.mod
	$(FIDASIM_COMPILER) $(CFLAGS) -c -I$(NETCDF_INCLUDE) $<

fidasim_linalg.mod fidasim_linalg.o: fidasim_linalg.f90
	$(FIDASIM_COMPILER) $(CFLAGS) -c $<

clean:
	-rm *.mod *.o fidasim

clean_debug:
	-rm *.mod *.o fidasim_debug
