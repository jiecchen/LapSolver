% script initjava.m
%
% (re-)initializes java classes
% called by init function

clear java
javaaddpath target/LapSolver.jar
javaaddpath deps/matlabcontrol-4.1.0.jar
import lapsolver.*;
import lapsolver.algorithms.*;
import lapsolver.lsst.*;
import lapsolver.generators.*;
import lapsolver.util.*;
import lapsolver.solvers.*;
import lapsolver.solvers.kelner.*;
