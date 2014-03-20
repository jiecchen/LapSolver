function x = lapAmdSolver(la,b)
% function x = lapAmdSolver(la,b)
%
% wrapping of lapWrapSolver around amdSolver
%
% Copyright Daniel Spielman, 2014, Yale University.
% This file of the LapSolve package.
% See source file for license information


default('b',[]);

if (isempty(b)),
    f = lapWrapSolver('amdSolver',la);
    x = f;
    return
else
  x = lapWrapSolver('amdSolver',la,b);
end

    
