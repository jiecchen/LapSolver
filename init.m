% script init.m
%
% Created by Dan Spielman on Jun 26, 2007
%
% sets up the correct paths

root = getroot;
addpath([root, 'matlab']);
addpath([root, 'matlab/solvers']);
addpath([root, 'matlab/generators']);
addpath([root, 'matlab/util']);
addpath([root, 'matlab/sparsify']);
addpath([root, 'matlab/cut']);
addpath([root, 'matlab/tests']);

clear java
javaaddpath build/LapSolver.jar
import lapsolver.*;
import lapsolver.algorithms.*;
import lapsolver.lsst.*;
