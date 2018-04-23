# -*- Makefile -*-

FC        = gfortran     # F90 compiler
OPTIM     = -g -C     #-O3       # Optimization flags; set to -g -C for useful debug information


# It should hopefully not be necessary to edit anything below
FITSDIR = /home/erik/Nedlastinger/cfitsio/lib      # Directory containing libcfitsio.a
LAPACK  = -L/usr/lib -llapack -lblas
HEALPIX = -L/home/erik/Nedlastinger/Healpix_3.31/lib_local
HEALINC =  /home/erik/Nedlastinger/Healpix_3.31/include_local
OUTPUT  = cmbspec
FFLAGS  = $(HEALPIX) $(LAPACK)

# List of source files to be compiled
OBJS    = math_tools.o spline_1D_mod.o rk_mod.o bs_mod.o ode_solver.o \
	  params.o time_mod.o cmbspec.o rec_mod.o evolution_mod.o


# Linking stage
cmbspec: $(OBJS)
	$(FC) $(FFLAGS) $(OPTIM) -o $(OUTPUT) $(OBJS) $(LAPACK) #-L$(HEALPIX)


# Dependencies
cmbspec.o       : time_mod.o rec_mod.o evolution_mod.o
time_mod.o      : params.o spline_1D_mod.o ode_solver.o
evolution_mod.o : time_mod.o ode_solver.o rec_mod.o
 

# Compilation of source files
%.o : %.f90
	$(FC) $(OPTIM) -I$(HEALINC) -c $<

# Clean-up command (write "make clean")
.PHONY: clean
clean:
	rm *.mod *.o *~ cmbspec
