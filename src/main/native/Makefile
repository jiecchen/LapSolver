CC=gcc
PROF_FLAGS=-fprofile-generate
#PROF_FLAGS=-fprofile-use
TARGET=../resources/lib
LIBRARY_NAME=liblapsolver
CFLAGS=-O3 -march=native -mtune=native -fopenmp -flto -std=c99 -fPIC
INCLUDES=-I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux -Iinclude -I/System/Library/Frameworks/JavaVM.framework/Versions/A/Headers
CLASSES=lapsolver.generators.Grid2
SOURCES=$(addsuffix .c,$(addprefix src/,$(subst .,_,$(CLASSES))))
OBJECTS=$(SOURCES:.c=.o)

LIB_NAME_OPT=-soname

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	DYNAMIC=-dynamiclib
	LIBRARY=$(addsuffix .jnilib,$(LIBRARY_NAME))
	CFLAGS+=-mno-avx
	INCLUDES=-I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/darwin -Iinclude
endif
ifeq ($(UNAME_S), Linux)
	DYNAMIC=-shared -Wl,-soname,liblapsolver
	LIBRARY=$(addsuffix .so,$(LIBRARY_NAME))
	INCLUDES=-I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux -Iinclude
endif

all: $(SOURCES) $(LIBRARY) tests

$(LIBRARY): $(OBJECTS)
	$(CC) $(CFLAGS) $(PROF_FLAGS) $(DYNAMIC) $^ -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(PROF_FLAGS) $(INCLUDES) -c $< -o $@

tests: grid2benchmark
grid2benchmark: tests/grid2benchmark.o $(OBJECTS)
	$(CC) $(CFLAGS) $(PROF_FLAGS) $^ -o $@

install: all
	mv $(LIBRARY) $(TARGET)/

.PHONY: clean

clean:
	rm -f src/*.o tests/*.o
	rm -f $(TARGET)/$(LIBRARY)

