#!/bin/bash
cd $1
# Set up library paths
source /opt/intel/composerxe/bin/compilervars.sh intel64
export DYLD_LIBRARY_PATH=../../:$DYLD_LIBRARY_PATH
export LD_LIBRARY_PATH=../../:$LD_LIBRARY_PATH
# Munch the arguments
PROG=$1
shift 1
# Run the test
if [ "$TESTMODE" = "auto" ]; then
	./$PROG $@ < test.in | cmp - test.out
else
	cp ../Makefile.test Makefile
	make 2>&1 > /dev/null
	$VALGRIND ./$PROG $@ < test.in
	[[ $KEEP = "yes" ]] || make clean 2>&1 > /dev/null
fi
