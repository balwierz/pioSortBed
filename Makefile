VERSION  = 3.0.12
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

# Build a self-contained Linux x86_64 release binary with libtbb statically
# linked. Requires building oneTBB locally first because Debian/Ubuntu only
# ship libtbb as a shared object (libtbb.so), not as an archive (libtbb.a).
#
# Set TBB_LIB to your locally-built libtbb.a, e.g.:
#   make release-binary TBB_LIB=/path/to/oneTBB/build/.../libtbb.a
# To build oneTBB statically:
#   git clone --branch v2022.3.0 https://github.com/uxlfoundation/oneTBB.git
#   cd oneTBB && mkdir build && cd build
#   cmake -DTBB_TEST=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release \
#         -DCMAKE_CXX_FLAGS='-fPIC -Wno-error -Wno-stringop-overflow' .. && make -j
TBB_LIB ?=
release-binary: src/pioSortBed.o
	@if [ -z "$(TBB_LIB)" ] || [ ! -f "$(TBB_LIB)" ]; then \
		echo "Error: set TBB_LIB to a static libtbb.a"; \
		echo "  e.g. make release-binary TBB_LIB=/path/to/libtbb.a"; \
		exit 1; \
	fi
	$(CC) -o pioSortBed-linux-x86_64 src/pioSortBed.o $(CFLAGS) \
	    -static-libstdc++ -static-libgcc "$(TBB_LIB)" -lpthread

clean:
	rm -f src/pioSortBed.o pioSortBed pioSortBed-linux-x86_64

.PHONY: install test release-binary clean
