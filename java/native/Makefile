CC=gcc
#PROF_FLAGS=-fprofile-generate
PROF_FLAGS=-fprofile-use
CFLAGS= -ftree-vectorizer-verbose=1 -fopenmp -flto -O3 -mtune=native -march=native -std=c99 -fPIC -g
INCLUDES=-I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux
TARGET=../../build

all: liblapsolver.so test

liblapsolver.so: lapsolver_generators_Grid2.o
	$(CC) $(CFLAGS) $(PROF_FLAGS) $^ -shared -o $@ -Wl,-soname,liblapsolver

test: test.o lapsolver_generators_Grid2.o
	$(CC) $(CFLAGS) $(PROF_FLAGS) $^ -o $@

test.o: test.c
	$(CC) $(CFLAGS) $(PROF_FLAGS) $(INCLUDES) -c $^

lapsolver_generators_Grid2.o: lapsolver_generators_Grid2.c
	$(CC) $(CFLAGS) $(PROF_FLAGS) $(INCLUDES) -c $^

install: all
	cp liblapsolver.so $(TARGET)/

clean:
	rm *.so *.o test

