% script initjava.m
%
% (re-)initializes java classes
% called by init function

clear java
javaaddpath build/LapSolver.jar
import lapsolver.*;
import lapsolver.algorithms.*;
import lapsolver.lsst.*;
import lapsolver.generators.*;
import lapsolver.util.*;