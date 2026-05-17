VERSION  = 3.8.0
PREFIX  ?= /usr/local

CC      = g++
CFLAGS  = -Isrc -O3 -std=c++20 -DVERSION_STRING=\"$(VERSION)\"
# Parallel std::sort on libstdc++ is implemented over TBB.
# C++/GCC runtimes statically linked so stdio doesn't go through a PLT;
# libtbb stays dynamic (no libtbb.a in Debian).
# liblz4 / libzstd back the --external-merge codecs (compressed binary temp
# files); both ship in libfoo.{a,so} on every distro that has -dev packages
# (`apt install liblz4-dev libzstd-dev`). Dynamic-linked here; the release
# binary target rolls them into a static build.
LDFLAGS = -static-libstdc++ -static-libgcc -ltbb -llz4 -lzstd
DEPS    =

# Optional htslib-backed features. Three independent opt-in flags,
# each enables a feature that links against htslib:
#   WITH_HTSLIB=1 — integrated bgzip + tabix output (--bgzip / --tabix),
#                   producing a queryable .bed.gz + .tbi in one step.
#                   Works with a system htslib install (`apt install
#                   libhts-dev`); HTSLIB= optional.
#   WITH_BAM=1    — BAM input (`*.bam` → coord-sorted BAM via htslib's
#                   BGZF worker pool). Implies WITH_HTSLIB.
#   WITH_RANS=1   — htscodecs SIMD rANS codecs for --external-merge temp
#                   files (--merge-codec=rans0|rans1). Requires HTSLIB
#                   to point at an htslib *source tree* (the htscodecs/
#                   headers ship there, not in the system install).
#                   Implies WITH_HTSLIB.
#
# Default (no flags): binary is BED-only with zero htslib dependency.
WITH_HTSLIB ?=
WITH_BAM    ?=
WITH_RANS   ?=
HTSLIB      ?=
ifneq ($(WITH_BAM),)
  WITH_HTSLIB := 1
endif
ifneq ($(WITH_RANS),)
  WITH_HTSLIB := 1
endif
ifneq ($(WITH_HTSLIB),)
  CFLAGS += -DWITH_HTSLIB
  ifneq ($(HTSLIB),)
    # Source-tree htslib: use its headers + statically link libhts.a so the
    # htscodecs rANS internals are pulled in (they're in libhts.a but not
    # exported by libhts.so).
    CFLAGS  += -I$(HTSLIB)
    LDFLAGS += $(HTSLIB)/libhts.a -lz -lpthread -ldeflate
  else
    # System htslib install: -lhts is enough for bgzip + tabix + BAM.
    LDFLAGS += -lhts -lz -lpthread -ldeflate
  endif
  ifneq ($(WITH_BAM),)
    CFLAGS += -DWITH_BAM
  endif
  ifneq ($(WITH_RANS),)
    ifeq ($(HTSLIB),)
      $(error WITH_RANS=1 needs HTSLIB=/path/to/htslib-source-tree (for the htscodecs/ headers))
    endif
    CFLAGS += -DWITH_RANS -I$(HTSLIB)/htscodecs
  endif
endif

# Optional LociSSD output (Parquet writer per FORMAT_SPEC.md). Build with:
#   make WITH_LOCISS=1
# Requires libarrow-dev + libparquet-dev (Debian/Ubuntu) or equivalent.
# Adds the --lociss-output FILE CLI option that emits a v2 LociSSD
# Parquet file instead of (or alongside) BED text. Default `make` (no
# WITH_LOCISS) has zero Arrow/Parquet dependency.
WITH_LOCISS ?=
ifneq ($(WITH_LOCISS),)
  CFLAGS  += -DWITH_LOCISS
  LDFLAGS += -lparquet -larrow
endif

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
	    -static-libstdc++ -static-libgcc "$(TBB_LIB)" \
	    -Wl,-Bstatic -llz4 -lzstd -Wl,-Bdynamic -lpthread

clean:
	rm -f src/pioSortBed.o pioSortBed pioSortBed-linux-x86_64

.PHONY: install test release-binary clean
