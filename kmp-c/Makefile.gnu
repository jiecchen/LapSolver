SHELL := /bin/bash
CC:=g++-4.9
BUILD:=build
SRC:=src
TARGET:=kmp
MKLROOT=/opt/intel/composerxe/mkl
MKL_CFLAGS = -flto -fopenmp -I$(MKLROOT)/include
MKL_LFLAGS = -flto -fopenmp
ifeq "$(LINKTYPE)" "static"
	MKL_LFLAGS += $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_intel_thread.a
	MKL_CFLAGS += -static
	MKL_LFLAGS += -static
else
	MKL_LFLAGS += -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread
endif
MKL_LFLAGS += -lpthread -lm

CFLAGS:=-fcilkplus -lcilkrts -O3 -march=native -mtune=native -std=c++11 -I$(SRC)
LFLAGS:=$(CFLAGS)

SOURCES := $(shell find $(SRC) -name *.cpp)
OBJECTS := $(addprefix $(BUILD)/,$(notdir $(SOURCES:.cpp=.o)))

.PHONY: all clean
all: $(BUILD) $(TARGET)
clean:
	rm -rf $(BUILD)
	rm -f $(TARGET)

$(BUILD):
	mkdir -p $@

$(TARGET): $(OBJECTS)
	$(CC) $(LFLAGS) $^ -o $@ $(MKL_LFLAGS)

define sourcerule
$(patsubst %.cpp, $(BUILD)/%.o, $(notdir $(1))): $(1)
	$$(CC) $$(CFLAGS) $$(MKL_CFLAGS) -c $$< -o $$@
endef

$(foreach src, $(SOURCES), $(eval $(call sourcerule, $(src))))
