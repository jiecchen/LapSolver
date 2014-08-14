CC:=icc
BUILD:=build
SRC:=src
TARGET:=libkmp.dylib
MKLROOT=/opt/intel/composerxe/mkl
MKL_CFLAGS = -ipo -openmp -mkl=parallel
MKL_LFLAGS = -ipo -openmp -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm

CFLAGS:=-O3 -xHost -std=c++11 -fPIC -I$(SRC)
LFLAGS:=$(CFLAGS)

SHARED_FLAG:=-dynamiclib

SOURCES := $(shell find $(SRC) -name *.cpp)
OBJECTS := $(addprefix $(BUILD)/,$(notdir $(SOURCES:.cpp=.o)))

.PHONY: all clean test
all: $(BUILD) $(TARGET)
clean:
	rm -rf $(BUILD)
	rm -f $(TARGET)
	rm -f test/$(TARGET)

test: $(BUILD) $(TARGET)
	@cp $(TARGET) test/$(TARGET)
	@echo
	@echo =======================
	@echo "        TESTING"
	@echo =======================
	@echo
	@make -C test
	@echo

$(BUILD):
	mkdir -p $@

$(TARGET): $(OBJECTS)
	$(CC) $(LFLAGS) $^ $(SHARED_FLAG) -o $@ $(MKL_LFLAGS)

define sourcerule
$(patsubst %.cpp, $(BUILD)/%.o, $(notdir $(1))): $(1)
	$$(CC) $$(CFLAGS) $$(MKL_CFLAGS) -c $$< -o $$@
endef

$(foreach src, $(SOURCES), $(eval $(call sourcerule, $(src))))
