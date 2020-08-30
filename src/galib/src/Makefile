FC      = gfortran
FCFLAGS = -O2 -Wall -Wno-unused-dummy-argument -Wno-unused-function -Wno-uninitialized -Wno-conversion -static -ffree-line-length-0 -std=gnu -fbacktrace -fbounds-check
# FCFLAGS = -O2 -Wall -Wno-unused-dummy-argument -Wno-unused-function# -static
# FCFLAGS = -O2 -ffpe-trap=invalid,zero,overflow -fbacktrace -g
# FCFLAGS = -fopenmp -static-libgfortran -static-libgcc

DEST_LIB = ../lib
DEST_INC = ../include

TARLIBD = libfec.dll
TARLIBS = libfec.a
OBJS    = moeadc.o moead.o moead_init.o moeadi.o moead_unit.o tnsdm.o \
          nsga2c.o nsga2.o nsga2_init.o nsga2i.o soga.o soga_init.o sogai.o\
          ga_unit.o ga_unit_selection.o ga_unit_crossover.o ga_unit_mutation.o \
          basic_optimizer.o individual.o kriging.o problem.o matrix.o util.o interface.o

LIBPATH = .
LIBNAME = fec
LIBOPT  = $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(LIBNAME))
INCOPT  = $(addprefix -I,$(LIBPATH))
LIBBLAS = -llapack -lrefblas
INCLUDE = ../include


################################################################################

.SUFFIXES: .f08
.PHONY: clean


all: st
dy: dynamic
st: static

dynamic: $(TARLIBD)
static: $(TARLIBS)

$(TARLIBD): $(OBJS)
	$(FC) $(FCFLAGS) -shared -o $@ $^ $(LIBOPT) $(LIBBLAS)
	mkdir -p $(DEST_LIB) $(DEST_INC)
	cp $@ $(DEST_LIB)
	cp *.a $(DEST_LIB)
	cp *.mod *.smod $(DEST_INC)

$(TARLIBS): $(OBJS)
	ar cr $@ $^
	mkdir -p $(DEST_LIB) $(DEST_INC)
	# cp $@ $(DEST_LIB)
	cp *.a $(DEST_LIB)
	cp *.mod *.smod $(DEST_INC)

%.o: %.f08
	$(FC) $(FCFLAGS) -c -o $@ $<

%.mod: %.f08 %.o
	@:

%.smod: %.f08 %.o
	@:

clean:
	rm -f *.o *.mod *.smod $(DEST_LIB)/* $(DEST_INC)/*


wrap: wrapper_xover.f08
	$(FC) $(FCFLAGS) -shared -o xover.dll $< $(LIBOPT) $(INCOPT)


################################################################################

###  common ###
matrix.o: util.o

### core ###
problem.o: interface.o individual.o
kriging.o: problem.o matrix.o util.o
ga_unit.o: interface.o util.o
ga_unit_selection: ga_unit.o
ga_unit_crossover.o: ga_unit.o
ga_unit_mutation.o: ga_unit.o
moead_unit.o: interface.o util.o
individual.o: interface.o util.o
basic_optimizer.o: interface.o individual.o problem.o

### driver ###
soga.o: interface.o util.o individual.o ga_unit.o problem.o basic_optimizer.o
soga_init.o: soga.o
sogai.o: soga.o
nsga2.o: interface.o util.o individual.o ga_unit.o problem.o soga.o
nsga2_init.o: nsga2.o
nsga2i.o: nsga2.o
nsga2c.o: interface.o util.o individual.o ga_unit.o problem.o soga.o nsga2.o
nsga2c_dp.o: interface.o util.o individual.o ga_unit.o problem.o soga.o nsga2.o nsga2c.o
tnsdm.o: interface.o util.o individual.o ga_unit.o problem.o soga.o nsga2c.o
tnsdm_1.o: tnsdm.o
moead.o: interface.o util.o individual.o moead_unit.o ga_unit.o problem.o basic_optimizer.o
moead_init.o: moead.o
moeadc.o: interface.o util.o individual.o moead_unit.o ga_unit.o problem.o moead.o
moeadi.o: moead.o

# ga_unit_selection.o ga_unit_crossover.o ga_unit_mutation.o

################################################################################

### TEST: Kriging ###

TARKRIG = krig.exe
OBJKRIG = test_kriging.o

krig: $(TARKRIG)

$(TARKRIG): $(OBJKRIG)
	$(FC) $(FCFLAGS) -o $@ $^ $(LIBOPT) $(LIBBLAS)

test_kriging.o: libfec.a


### TEST: lib ###

TARTEST = tes.exe
OBJTEST = test_tnsdm.o
# OBJTEST = test_nsga2.o

test: $(TARTEST)

$(TARTEST): $(OBJTEST)
	$(FC) $(FCFLAGS) -o $@ $^ $(LIBOPT)

test_tnsdm.o: libfec.a
