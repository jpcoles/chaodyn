CC=gcc
OPENMP=-fopenmp
CFLAGS=-Wall -g -O3 -std=c99 $(OPENMP)
CFLAGS+=-ftree-vectorizer-verbose=3 -ftree-vectorize
CFLAGS+=-fno-omit-frame-pointer -floop-optimize2 -funroll-loops -fprefetch-loop-arrays
CFLAGS+=-fstrict-aliasing -mpreferred-stack-boundary=4 

LDFLAGS=-lpng -lm $(OPENMP)

all: mss3domp

tdyn: io.o tdyn.o
	$(CC) $(LDFLAGS) $^ -o $@ 

sim2png: io.o sim2png.o
	$(CC) $(LDFLAGS) $^ -o $@ 

mss3domp: io.o mss3domp.o
	$(CC) $(LDFLAGS) $^ -o $@ 

simdiff: simdiff.o
	$(CC) $(LDFLAGS) $^ -o $@ 

mss: mss.o
	$(CC) $(LDFLAGS) $^ -o $@ 

mss3d: mss3d.o
	$(CC) $(LDFLAGS) $^ -o $@ 

%.o: %.c
	$(CC) $(CFLAGS) -c $<
