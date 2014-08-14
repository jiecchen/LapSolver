TESTS:=$(patsubst %/,%,$(shell ls -d */))

.PHONY: all $(TESTS)

all: $(TESTS)

$(TESTS):
	@cp Makefile.test $@/Makefile
	@cd $@ && make 2>&1 > /dev/null && cd ..
	@./runtest.sh $@ ; [[ $$? -eq 0 ]] && echo "$@: passed" || echo "$@: failed"
	@cd $@ && make clean 2>&1 > /dev/null && cd ..
