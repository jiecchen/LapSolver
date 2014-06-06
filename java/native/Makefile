CC=gcc
PROF_FLAGS=-fprofile-use #-fprofile-generate
CFLAGS=-O3 -mtune=native -march=native -std=c99 -fPIC
INCLUDES=-I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux
TARGET=../../build

liblapsolver.so: lapsolver_generators_Grid2.o
	$(CC) $(PROF_FLAGS) $^ -shared -o $@ -Wl,-soname,liblapsolver

lapsolver_generators_Grid2.o: lapsolver_generators_Grid2.c
	$(CC) $(CFLAGS) $(PROF_FLAGS) $(INCLUDES) -c $^

install: liblapsolver.so
	cp liblapsolver.so $(TARGET)/

clean:
	rm *.so *.o
