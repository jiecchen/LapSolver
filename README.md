# LapSolver
An ultra-fast solver for Laplacian systems of linear equations.

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

## Instructions

How to prepare LapSolver on OS X:

    1) Install JDK7+ from Oracle's website
    2) Set $JAVA_HOME to the JDK8 directory
        You can do this by adding:

        > export JAVA_HOME=`/usr/libexec/java_home`

        to ~/.profile
    3) Using Homebrew or Macports, install Maven
    4) ** optional ** Using Homebrew or Macports, install GCC 4.8+
    5) From the top level folder, run:

        > mvn install

  For Matlab:

    1) Requires R2014a, or an earlier version equipped with JDK7+

  For Python (Anaconda):

    1) Install jpype — this requires modifying the included setup.py.
       The default setupMacOSX function should begin:

               def setupMacOSX(self):
                   self.javaHome = '{{VALUE OF JAVA_HOME}}'
                   self.jdkInclude = "darwin"

    2) python setup.py install
    3) import lapsolver!
