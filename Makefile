CC=g++
CFLAGS=-I.
# static can be removed from the options, but it is safer to keep it and do not recompile often.
LDFLAGS=-lboost_program_options -static
DEPS = 

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

pioSortBed: pioSortBed.o
	$(CC) -o $@ pioSortBed.o $(CFLAGS) $(LDFLAGS)
