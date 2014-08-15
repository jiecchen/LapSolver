cd $1
# Set up library paths
source /opt/intel/composerxe/bin/compilervars.sh intel64
export DYLD_LIBRARY_PATH=../:$DYLD_LIBRARY_PATH
export LD_LIBRARY_PATH=../:$LD_LIBRARY_PATH
# Run the test
if [ "$TESTMODE" = "auto" ]; then
	./$1 input | cmp - output
else
	cp ../Makefile.test Makefile
	make 2>&1 > /dev/null
	./$1 input
	make clean 2>&1 > /dev/null
fi
