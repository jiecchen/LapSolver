CC:=icc
TARGET:=kmp
MKLROOT=/opt/intel/composerxe/mkl
MKL_CFLAGS = -ipo -openmp -mkl=parallel
MKL_LFLAGS = -ipo -openmp
ifeq "$(LINKTYPE)" "static"
	MKL_LFLAGS += $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_intel_thread.a
	MKL_CFLAGS += -static-intel
	MKL_LFLAGS += -static-intel
else
	MKL_LFLAGS += -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread
endif
MKL_LFLAGS += -lpthread -lm

CFLAGS:=-O3 -std=c++11 -mtune=native -march=native -I.
LFLAGS:=$(CFLAGS)

SOURCES = main.cpp structures/edgelist.cpp
OBJECTS=$(SOURCES:.cpp=.o)

.PHONY: all clean
all: $(TARGET)
clean:
	rm -f $(TARGET)
	rm -f $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(LFLAGS) $^ -o $@ $(MKL_LFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(MKL_CFLAGS) -c $< -o $@


