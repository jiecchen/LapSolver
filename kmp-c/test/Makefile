SHELL := /bin/bash
TESTS:=$(patsubst %/,%,$(shell ls -d */))
BENCH:=$(addsuffix -bench, $(TESTS))

.PHONY: all bench $(TESTS) $(BENCH)

all: $(TESTS)
bench: $(BENCH)

$(TESTS):
	@cp Makefile.test $@/Makefile
	@cd $@ && make 2>&1 > /dev/null && cd ..
	@TESTMODE=auto ./runTest.sh $@ ; [[ $$? -eq 0 ]] && echo "$@: passed" || echo "$@: failed"
	@cd $@ && make clean 2>&1 > /dev/null && cd ..

$(BENCH):
	$(eval TEST := $(patsubst %-bench, %, $@))
	@cp Makefile.test $(TEST)/Makefile
	@cd $(TEST) && make 2>&1 > /dev/null && cd ..
	@TESTMODE=auto ./runTest.sh $(TEST) -v -n 50
	@cd $(TEST) && make clean 2>&1 > /dev/null && cd ..
