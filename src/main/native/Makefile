CC=gcc
PROF_FLAGS=-fprofile-generate
#PROF_FLAGS=-fprofile-use
TARGET=../resources/lib
LIBRARY=liblapsolver.so
CFLAGS=-O3 -mtune=native -march=native -fopenmp -flto -std=c99 -fPIC
INCLUDES=-I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux -Iinclude
CLASSES=lapsolver.generators.Grid2
SOURCES=$(addsuffix .c,$(addprefix src/,$(subst .,_,$(CLASSES))))
OBJECTS=$(SOURCES:.c=.o)

all: $(SOURCES) $(LIBRARY) tests

$(LIBRARY): $(OBJECTS)
	$(CC) $(CFLAGS) $(PROF_FLAGS) $^ -shared -o $@ -Wl,-soname,liblapsolver

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

