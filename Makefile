CC=gcc
OPENMP=#-fopenmp
INCDIR=-Itipsylib/trunk
INCDIR+=-I/usr/X11/include
LIBDIR=-Ltipsylib/trunk
LIBDIR+=-L/usr/X11/lib
CFLAGS=-Wall -g -O3 -std=c99 $(OPENMP) $(INCDIR)
CFLAGS+=-ftree-vectorizer-verbose=3 -ftree-vectorize
CFLAGS+=-fno-omit-frame-pointer #-floop-optimize2 
CFLAGS+=-funroll-loops -fprefetch-loop-arrays
CFLAGS+=-fstrict-aliasing -mpreferred-stack-boundary=4 
CFLAGS+=-fnested-functions

LDFLAGS=$(LIBDIR) -lpng -lm $(OPENMP) -ltipsy

default: chaodyn1D

all: chaodyn tdyn sim2png simdiff chaodyn1D simcat simprop

chaodyn: io.o io_tipsy.o ic.o chaodyn.o
	$(CC) $(LDFLAGS) $^ -o $@ 

simcat: io.o simcat.o
	$(CC) $(LDFLAGS) $^ -o $@ 

tdyn: io.o tdyn.o
	$(CC) $(LDFLAGS) $^ -o $@ 

sim2png: io.o sim2png.o
	$(CC) $(LDFLAGS) $^ -o $@ 

simprop: io.o simprop.o
	$(CC) $(LDFLAGS) $^ -o $@ 

simdiff: simdiff.o
	$(CC) $(LDFLAGS) $^ -o $@ 

chaodyn1D: chaodyn1D.o
	$(CC) $(LDFLAGS) $^ -o $@ 

%.o: %.c
	$(CC) $(CFLAGS) -c $<


clean:
	-rm *.o chaodyn tdyn sim2png simdiff chaodyn1D
