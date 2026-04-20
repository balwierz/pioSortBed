VERSION  = 2.0.0
PREFIX  ?= /usr/local

CC      = g++
CFLAGS  = -Isrc -O3 -std=c++11 -fopenmp -DVERSION_STRING=\"$(VERSION)\"
# static can be removed from the options, but it is safer to keep it and do not recompile often.
LDFLAGS = -lgomp
DEPS    =

src/pioSortBed.o: src/pioSortBed.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

pioSortBed: src/pioSortBed.o
	$(CC) -o $@ src/pioSortBed.o $(CFLAGS) $(LDFLAGS)

install: pioSortBed
	install -m 755 pioSortBed $(PREFIX)/bin/pioSortBed

clean:
	rm -f src/pioSortBed.o pioSortBed

.PHONY: install clean
