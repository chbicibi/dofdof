FC      = gfortran
FCFLAGS = -O2 -Wall -Wno-unused-dummy-argument -Wno-unused-function -static
# FCFLAGS = -O2 -ffpe-trap=invalid,zero,overflow -fbacktrace -g
# FCFLAGS = -fopenmp -static-libgfortran -static-libgcc

DEST_LIB = ../lib
DEST_INC = ../include

TARLIBD = libkriging.dll
TARLIBS = libkriging.a
OBJS    = kriging.o matrix.o util.o

LIBPATH = .
LIBNAME = kriging
INCPATH = ../include
LIBOPT  = $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(LIBNAME))
INCOPT  = $(addprefix -I,$(INCPATH))
LIBBLAS = -llapack -lrefblas


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
	# cp *.a $(DEST_LIB)
	cp *.mod $(DEST_INC)

$(TARLIBS): $(OBJS)
	ar cr $@ $^
	mkdir -p $(DEST_LIB) $(DEST_INC)
	cp $@ $(DEST_LIB)
	# cp *.a $(DEST_LIB)
	cp *.mod $(DEST_INC)

%.o: %.f08
	$(FC) $(FCFLAGS) -c -o $@ $<

%.mod: %.f08 %.o
	@:

%.smod: %.f08 %.o
	@:

clean:
	rm -f *.o *.mod $(DEST_LIB)/* $(DEST_INC)/*


################################################################################

matrix.o: util.o

kriging.o: matrix.o util.o

util.o:
	$(FC) $(FCFLAGS) -c -o $@ ../util.f08

################################################################################
