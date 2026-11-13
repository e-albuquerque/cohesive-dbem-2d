################################################################################
# MAKEFILE - Just edit if you know what you are doing!! 
################################################################################

# Compiler
FC = gfortran 

# Flags for debugging or for maximum performance (comment as necessary)
FLAGS1 = -O3 -fopenmp -fimplicit-none
# FLAGS1 = -g -O0 -fbounds-check -fopenmp -fimplicit-none # Debug mode 
FCFLAGS = $(FLAGS1) 

# Libraries needed for linking
LDFLAGS = -llapack 
#LDFLAGS = liblapack.so.3.4.2
#LIBS = -L/usr/local/lib/gnu -llapack -lblas #CENAPAD!!
#LIBS = -L/usr/lib/x86_64-linux-gnu -llapack
LIBS = -llapack -lblas

# List of executables to be built 
PROGRAMS = MAIN_COHESIVE_DBEM_2D

# "make" builds all
all: $(PROGRAMS)	

MAIN_COHESIVE_DBEM_2D: MODULE_DATA.o MODULE_INPUT_DATA.o MODULE_DATA_TREATMENT.o MODULE_ANALYSIS.o MODULE_POSTPROCESSING.o MODULE_SOLVER.o MODULE_APPROXIMATION_FUNCTIONS.o MODULE_GEOMETRIC_MAPPINGS.o MODULE_NUMERICAL_INTEGRATION.o

MAIN_COHESIVE_DBEM_2D.o: MODULE_DATA.o MODULE_INPUT_DATA.o MODULE_DATA_TREATMENT.o MODULE_ANALYSIS.o MODULE_POSTPROCESSING.o MODULE_SOLVER.o MODULE_APPROXIMATION_FUNCTIONS.o MODULE_GEOMETRIC_MAPPINGS.o MODULE_NUMERICAL_INTEGRATION.o

MODULE_APPROXIMATION_FUNCTIONS.o: MODULE_DATA.o

MODULE_GEOMETRIC_MAPPINGS.o: MODULE_DATA.o

MODULE_NUMERICAL_INTEGRATION.o: MODULE_DATA.o MODULE_GEOMETRIC_MAPPINGS.o MODULE_APPROXIMATION_FUNCTIONS.o

MODULE_INPUT_DATA.o: MODULE_DATA.o MODULE_NUMERICAL_INTEGRATION.o

MODULE_DATA_TREATMENT.o: MODULE_DATA.o MODULE_GEOMETRIC_MAPPINGS.o

MODULE_SOLVER.o: MODULE_DATA.o MODULE_GEOMETRIC_MAPPINGS.o

MODULE_ANALYSIS.o: MODULE_DATA.o MODULE_NUMERICAL_INTEGRATION.o MODULE_SOLVER.o

MODULE_POSTPROCESSING.o: MODULE_DATA.o MODULE_GEOMETRIC_MAPPINGS.o MODULE_APPROXIMATION_FUNCTIONS.o

# $@: whatever the target is (anything to the left of the ":")
# $^: whatever the target dependes on (anything to the right of the ":")
# $<: used to list only the first pre-requisite (source file)

# Building prog from prog.o; $^ (GNU extension) is
%: %.o
	@echo 'Building target: $@'
	@echo 'Invoking: GNU Fortran Linker'
	$(FC) $(FCFLAGS) -o $@ $^ $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Building prog.o from prog.f90 
%.o: %.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	$(FC) $(FCFLAGS) -c $< $(LDFLAGS)
	@echo 'Finished building: $<'
	@echo ' '

# Utility targets
.PHONY: clean veryclean

clean:
	-rm -f *.o *.mod *.MOD

veryclean: clean
	-rm -f *~ $(PROGRAMS) $(PROGRAMS).exe
