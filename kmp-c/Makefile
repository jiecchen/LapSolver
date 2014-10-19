SHELL:=/bin/bash
CC:=icc
BUILD:=build

TARGET:=libkmp.dylib
ifeq ($(shell uname -s), Linux)
	TARGET:=libkmp.so	
endif

.PHONY: all test bench clean

all: $(TARGET)

test: all
	@echo
	@echo =======================
	@echo "        TESTING"
	@echo =======================
	@echo
	@make -s -C test
	@echo

bench: all
	@echo
	@echo =======================
	@echo "       BENCHMARK"
	@echo =======================
	@echo
	@make -s -C test bench
	@echo

clean:
	rm -rf $(TARGET) $(BUILD)

$(BUILD):
	mkdir -p $@

$(TARGET): $(BUILD)
	cd $(BUILD) && CC=${CC} CXX=${CC} cmake ..
	cd $(BUILD) && make
	cp $(BUILD)/$(TARGET) .
