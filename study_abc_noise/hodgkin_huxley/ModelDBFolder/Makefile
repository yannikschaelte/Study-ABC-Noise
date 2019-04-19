# Modified from: http://www.webalice.it/o.drofa/davide/makefile-fortran/makefile-fortran.html
# Create Executable by typing make in command line

# ======================================================================
# Start with the declarations
# ======================================================================

# The compiler
FC = gfortran

# flags for debugging or for maximum performance, comment as necessary
#FCFLAGS = -g -fbounds-check
#FCFLAGS = -O2

# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/include

# libraries needed for linking
LDFLAGS = -L/usr/local/epd/lib -I/usr/local/epd/include -llapack -lblas

# List of executables to be built within the package
PROGRAMS = HH_run

# "make" builds all
all: $(PROGRAMS)

# ======================================================================
# Rules for programs, modify as needed
# ======================================================================

HH_run.o: ParameterModule.o MT19937.o HH_master.o 
HH_run: ParameterModule.o MT19937.o HH_master.o 

# ======================================================================
# General rules, should not require modification
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f95
	$(FC) $(FCFLAGS) -c $<

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)

