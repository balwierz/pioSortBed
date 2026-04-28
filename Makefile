VERSION  = 2.1.0
PREFIX  ?= /usr/local

CC      = g++
CFLAGS  = -Isrc -O3 -std=c++17 -DVERSION_STRING=\"$(VERSION)\"
# Parallel std::sort on libstdc++ is implemented over TBB.
# C++/GCC runtimes statically linked so stdio doesn't go through a PLT;
# libtbb stays dynamic (no libtbb.a in Debian).
LDFLAGS = -static-libstdc++ -static-libgcc -ltbb
DEPS    =

src/pioSortBed.o: src/pioSortBed.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

pioSortBed: src/pioSortBed.o
	$(CC) -o $@ src/pioSortBed.o $(CFLAGS) $(LDFLAGS)

install: pioSortBed
	install -m 755 pioSortBed $(PREFIX)/bin/pioSortBed

test: pioSortBed
	bash test/test.sh

clean:
	rm -f src/pioSortBed.o pioSortBed

.PHONY: install test clean
