cd $1
# Set up library paths
source /opt/intel/composerxe/bin/compilervars.sh intel64
export DYLD_LIBRARY_PATH=../../:$DYLD_LIBRARY_PATH
export LD_LIBRARY_PATH=../../:$LD_LIBRARY_PATH
# Run the test
if [ "$TESTMODE" = "auto" ]; then
	cat test.in | ./$1 | cmp - test.out
else
	cp ../Makefile.test Makefile
	make 2>&1 > /dev/null
	cat test.in | ./$1
	make clean 2>&1 > /dev/null
fi
