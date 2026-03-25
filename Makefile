CC=g++
CFLAGS=-I. -O3 -std=c++11 -fopenmp
# static can be removed from the options, but it is safer to keep it and do not recompile often.
LDFLAGS=-lgomp -static
DEPS = 

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

pioSortBed: pioSortBed.o
	$(CC) -o $@ pioSortBed.o $(CFLAGS) $(LDFLAGS)
