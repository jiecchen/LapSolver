SHELL := /bin/bash
CC:=icc
BUILD:=build
SRC:=../../src
TARGET:=$(shell basename $$(pwd))
MKLROOT=/opt/intel/composerxe/mkl
MKL_CFLAGS = -ipo -openmp -mkl=parallel
MKL_LFLAGS = -ipo -openmp
MKL_LFLAGS += -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm

KMP_LFLAGS += -L../../ -lkmp

CFLAGS:=-O3 -g -xHost -std=c++11 -I$(SRC)
LFLAGS:=$(CFLAGS)

LIBS:=$(MKL_LFLAGS) $(KMP_LFLAGS)

SOURCES := test.cpp
OBJECTS := $(addprefix $(BUILD)/,$(notdir $(SOURCES:.cpp=.o)))

.PHONY: all clean
all: $(BUILD) $(TARGET)
clean:
	rm -rf $(BUILD)
	rm -f $(TARGET)
	rm -f Makefile

$(BUILD):
	mkdir -p $@

$(TARGET): $(OBJECTS)
	$(CC) $(LFLAGS) $^ -o $@ $(LIBS)

define sourcerule
$(patsubst %.cpp, $(BUILD)/%.o, $(notdir $(1))): $(1)
	$$(CC) $$(CFLAGS) $$(MKL_CFLAGS) -c $$< -o $$@
endef

$(foreach src, $(SOURCES), $(eval $(call sourcerule, $(src))))
