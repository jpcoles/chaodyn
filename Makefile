CC=gcc
OPENMP=-fopenmp
CFLAGS=-Wall -g -O3 -std=c99 $(OPENMP)
CFLAGS+=-ftree-vectorizer-verbose=3 -ftree-vectorize
CFLAGS+=-fno-omit-frame-pointer -floop-optimize2 -funroll-loops -fprefetch-loop-arrays
CFLAGS+=-fstrict-aliasing -mpreferred-stack-boundary=4 

LDFLAGS=-lpng -lm $(OPENMP)

default: chaodyn

all: chaodyn tdyn sim2png simdiff chaodyn1D

tdyn: io.o tdyn.o
	$(CC) $(LDFLAGS) $^ -o $@ 

sim2png: io.o sim2png.o
	$(CC) $(LDFLAGS) $^ -o $@ 

chaodyn: io.o ic.o chaodyn.o
	$(CC) $(LDFLAGS) $^ -o $@ 

simdiff: simdiff.o
	$(CC) $(LDFLAGS) $^ -o $@ 

chaodyn1D: chaodyn1D.o
	$(CC) $(LDFLAGS) $^ -o $@ 

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	-rm *.o chaodyn tdyn sim2png simdiff chaodyn1D
