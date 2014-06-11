# LapSolver
An ultra-fast solver for LaPlacian systems of linear equations.

## Dependencies

  * Maven 3
  * Java 7
  * gcc 4.8 or clang equivalent
  * Mac OS X 10.8 or a recent Linux installation
  * Matlab 2014a to use from Matlab
  * Python 2.7.5 to use Python bindings

## How to use

Run "mvn install" to:

  * Compile Java/C code (called via JNI)
  * Create a JAR file for use with Matlab and other projects
  * Generate a Python package with bindings (located in the python/ directory)  

To use from Matlab, run "init" from this directory.
