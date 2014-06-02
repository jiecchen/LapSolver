function x = voltSolve(aflow,solver)
% function x = voltSolve(aflow,solver)
% function x = voltSolve(laflow,solver)
%
% for solving flow problems produced by flowProblem,
% using irls on the voltages.
%
% the key is that it sets up potentials of node 1 to 1
% of node 2 to 0,
% and then solves the system using solver, like amdSolver or such
%
% can take aflow matrix, or a laplacian version

if (min(min(sign(triu(aflow,1)))) == 0)
    la = lap(aflow);
else
    la = aflow;
end

lah = la(3:end,3:end);
b = la(3:end,1);

xh = solver(lah,-b);

x = [1;0;xh];
