CC=gcc
#PROF_FLAGS=-fprofile-generate
#PROF_FLAGS=-fprofile-use
TARGET=../resources/lib
LIBRARY_NAME=liblapsolver
CFLAGS=-O3 -march=native -mtune=native -flto -std=c99 -fPIC
LFLAGS:=$(CFLAGS)
CLASSES=lapsolver.generators.Grid2 \
        lapsolver.generators.TriangleGrid2
SOURCES=src/generators.c $(addsuffix .c,$(addprefix src/,$(subst .,_,$(CLASSES))))
OBJECTS=$(SOURCES:.c=.o)

ISCLANG:=$(shell $(CC) --version 2>&1 | grep clang ; echo $$?)
ifeq "$(ISCLANG)" "1"
	CFLAGS += -fopenmp
	LFLAGS += -fopenmp
endif

PLATFORM := $(shell uname -s | tr A-Z a-z)
ifeq ($(PLATFORM),darwin)
	DYNAMIC=-dynamiclib
	LIBRARY=$(addsuffix .dylib,$(LIBRARY_NAME))
	CFLAGS+=-mno-avx
endif
ifeq ($(PLATFORM), linux)
	DYNAMIC=-shared -Wl,-soname,$(LIBRARY_NAME)
	LIBRARY=$(addsuffix .so,$(LIBRARY_NAME))
endif

INCLUDES=-I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/$(PLATFORM) -Iinclude

.PHONY: all install $(CLASSES) clean profile-build

all: $(CLASSES) $(SOURCES) $(LIBRARY) tests

profile-build:
	make clean
	PROF_FLAGS=-fprofile-generate make all
	./bench.sh grid2benchmark
	make clean
	PROF_FLAGS=-fprofile-use make all

$(CLASSES):
	cd ../java && javah -d ../native/include $@

$(LIBRARY): $(OBJECTS)
	$(CC) $(LFLAGS) $(PROF_FLAGS) $(DYNAMIC) $^ -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(PROF_FLAGS) $(INCLUDES) -c $< -o $@

tests: grid2benchmark
grid2benchmark: tests/grid2benchmark.o $(OBJECTS)
	$(CC) $(LFLAGS) $(PROF_FLAGS) $^ -o $@

install: all
	mkdir -p $(TARGET)/
	mv $(LIBRARY) $(TARGET)/

clean:
	rm -f src/*.o tests/*.o grid2benchmark
	rm -f $(TARGET)/$(LIBRARY)

