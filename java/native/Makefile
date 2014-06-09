CC=gcc
PROF_FLAGS=#-fprofile-generate #-fprofile-use
CFLAGS=-flto -O3 -mtune=native -march=native -std=c99 -fPIC
INCLUDES=-I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux
TARGET=../../build

all: liblapsolver.so test

liblapsolver.so: lapsolver_generators_Grid2.o
	$(CC) $(CFLAGS) $(PROF_FLAGS) $^ -shared -o $@ -Wl,-soname,liblapsolver

test: test.c
	$(CC) $(CFLAGS) $^ -o $@

lapsolver_generators_Grid2.o: lapsolver_generators_Grid2.c
	$(CC) $(CFLAGS) $(PROF_FLAGS) $(INCLUDES) -c $^

install: all
	cp liblapsolver.so $(TARGET)/

clean:
	rm *.so *.o test

