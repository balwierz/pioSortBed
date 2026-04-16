CC=g++
CFLAGS=-Isrc -O3 -std=c++11 -fopenmp
# static can be removed from the options, but it is safer to keep it and do not recompile often.
LDFLAGS=-lgomp -static
DEPS =

src/pioSortBed.o: src/pioSortBed.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

pioSortBed: src/pioSortBed.o
	$(CC) -o $@ src/pioSortBed.o $(CFLAGS) $(LDFLAGS)
