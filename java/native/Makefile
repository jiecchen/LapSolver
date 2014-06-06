CC=gcc
CFLAGS=-O3 -g -ftree-vectorizer-verbose=1 -march=native -std=c99 -fPIC
INCLUDES=-I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux
TARGET=../../build

liblapsolver.so: lapsolver_generators_Grid2.o
	$(CC) $^ -shared -o $@ -Wl,-soname,liblapsolver

lapsolver_generators_Grid2.o: lapsolver_generators_Grid2.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $^

install: liblapsolver.so
	cp liblapsolver.so $(TARGET)/

clean:
	rm *.so *.o
