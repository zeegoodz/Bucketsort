
DEBUG=2
CC=mpicc

#CFLAGS=-Wall -g
CFLAGS=-Wall -O -DDEBUG=$(DEBUG)

#linker flags for libraries
LDFLAGS=-lprand
PROGS=parallel-bucketsort


all: $(PROGS)

parallel-bucketsort:	parallel-bucketsort.o timing.o
	$(CC) $(CFLAGS) -o $@ parallel-bucketsort.o timing.o $(LDFLAGS)

clean:
	/bin/rm --force *.o a.out $(PROGS)
