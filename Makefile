
DEBUG=2
CC=gcc

#CFLAGS=-Wall -g
CFLAGS=-Wall -O -DDEBUG=$(DEBUG)

#linker flags for libraries
LDFLAGS=
PROGS=sequential-bucketsort


all: $(PROGS)

sequential-bucketsort:	sequential-bucketsort.o timing.o
	$(CC) $(CFLAGS) -o $@ sequential-bucketsort.o timing.o $(LDFLAGS)

clean:
	/bin/rm --force *.o a.out $(PROGS)
