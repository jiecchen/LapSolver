function [cutval, v, stats, mats] = irlsCut(aflow,iters,solver,v0)
% function [cutval, v, stats, mats] = irlsCut(aflow,iters,solver,v0)
%
% use irls, as explained by anup, to
% find voltages of the min cut between node 1 and node 2
%
% assuming a graph produced by flowProblem
%
% runs for iters runs,
% takes a solver (e.g. amdSolver)
% and can take an initial vector
%
% right now stats is just a list of the solve times
% mats is a list of the matrices that needed to be solved
%
% does not presently handle wts
% Daniel A. Spielman, Yale University, Oct 23, 2012
%
% compare to:
% addpath ~/progs/matlab_bgl
% flowval = max_flow(aflow,1,2)

n = length(aflow);

[ai,aj] = find(tril(aflow));
m = e2m([ai,aj],-1);
la = m*m';

ne = size(m,2);

defaultStr('solver','@(la,b)(iccSolver(la,b))');

w = ones(ne,1);

if exist('v0')
  vd = (v0(:))'*m;
  w = 1./sqrt(abs(vd).^2 + 1e-5);
end

nmats = 1;

for iter = 1:iters,
    lap = m*diag(sparse(w))*m';
    
    tic;
    mats{nmats} = lap;
    nmats = nmats+1;
    v = voltSolve(lap,solver);
    stats(iter) = toc;
    
    vd = (v)'*m;
    w = 1./sqrt(abs(vd).^2 + 1e-5);

    
    [wt] = minVecCut(aflow,v)
    
end

cutval = wt;
