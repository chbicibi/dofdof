FC      = gfortran
# FCFLAGS = -O2 -Wall -Wno-unused-dummy-argument -Wno-unused-function -static
# FCFLAGS = -O2 -Wall -Wno-unused-dummy-argument -Wno-unused-function -Wno-uninitialized -Wno-conversion -static -ffree-line-length-0 -std=gnu -fbacktrace -fbounds-check
FCFLAGS = -O2 -Wall -static -ffree-line-length-0 -fbacktrace -fbounds-check
PARA    = -fopenmp
# FCFLAGS = -O2 -ffpe-trap=invalid,zero,overflow -fbacktrace -g -pedantic -std=f2008
# FCFLAGS = -fopenmp -static-libgfortran -static-libgcc

DEST     = ../

TARGET  = dof_new.exe
OBJS    = dof_kriging.o dof_base.o dof_flight.o dof_new.o

INCLUDE = galib/include
LIBPATH = galib/lib
LIBNAME = fec
INCOPT  = $(addprefix -I,$(INCLUDE))
LIBOPT  = $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(LIBNAME))
LIBBLAS = -llapack -lrefblas


################################################################################

.SUFFIXES: .f08
.PHONY: clean


all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FCFLAGS) $(PARA) -o $@ $^ $(LIBOPT) $(LIBBLAS)
	cp $@ $(DEST)

%.o: %.f08
	$(FC) $(FCFLAGS) $(PARA) -c -o $@ $< $(INCOPT)

%.mod: %.f08 %.o
	@:

%.smod: %.f08 %.o
	@:

clean:
	rm -f *.o *.mod *.smod


################################################################################

###  common ###
dof_flight.o: dof_kriging.o dof_base.o
dof_new.o: dof_kriging.o dof_base.o dof_flight.o
