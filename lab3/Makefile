################################################################
## SIMPLE MAKEFILE TO BUILD WITH MPICH                        ##
################################################################
# oct-2010 - a.knirsch@fbi.h-da.de                             #
# April-2013 - ronald.moore@h-da.de  (added make test)         #
################################################################

CXX	= mpic++
CFLAGS	= -O1
PROGS	= matrix_sequential matrix_parallel matrix_generator
OBJS	= matrix_sequential.o matrix_parallel.o matrix_generator.o CMatrix.o
TESTFILES = matrix-A-10x10.txt matrix-B-10x10.txt

PLATFORM=$(shell uname)
ifeq ($(PLATFORM), Darwin)
CFLAGS	+= -I/opt/local/include/
LDFLAGS += -L/opt/local/lib/
else
LDFLAGS +=
endif

all: $(PROGS) $(TESTFILES)

#      "Variables", die hier benutzt werden:
#      Vgl. http://www.makelinux.net/make3/make3-CHP-2-SECT-2.html
# $@ = The filename representing the target.
# $< = The filename of the first prerequisite.
# $* = The stem of the filename of the target (i.e. without .o, .cpp...)
# $^ = The names of all the prerequisites, with spaces between them.

$(PROGS): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@  $@.o  CMatrix.o
#	$(CXX) $(LDFLAGS) -o matrix_generator  matrix_generator.o  CMatrix.o
#	$(CXX) $(LDFLAGS) -o matrix_sequential matrix_sequential.o CMatrix.o
#	$(CXX) $(LDFLAGS) -o matrix_parallel   matrix_parallel.o   CMatrix.o

%.o: %.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

$(TESTFILES): matrix_generator
	./matrix_generator 10 10 >$@

clean:
	rm -rf $(PROGS) *~ *.o ._*

distclean: clean
	rm -rf $(TESTFILES)

test:  $(PROGS) $(TESTFILES)
	echo "Calling matrix_sequential..."
	./matrix_sequential $(TESTFILES)
	echo "Calling matrix_parallel..."
	mpiexec ./matrix_parallel $(TESTFILES)
