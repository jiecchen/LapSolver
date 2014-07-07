function aflow = flowProblem(a,k)
% function aflow = flowProblem(a,k)
%
% takes as input a graph a (in adj format).
% creates two new nodes, 1 and 2,
% and the problem is to find the maximum flow between them.
%
% each of these nodes are connected to bfs balls containing k
% vertices
%
% the code tries to put them as far apart as possible
% but, there can be overlap
%

import lapsolver.*;
import lapsolver.util.*;

n = length(a);

[ai,aj,av] = find(tril(a));

start = ceil(n*rand(1));

g = Graph(ai-1,aj-1, av);
bfs = 1 + GraphUtils.bfsOrder(g,start-1);

s1 = bfs(1:k);

last = bfs(end);
bfs = 1 + GraphUtils.bfsOrder(g,last-1);

s2 = bfs(1:k);

aflow((3:n+2),(3:n+2)) = a;
aflow(1,s1+2) = 1;
aflow(s1+2,1) = 1;
aflow(2,s2+2) = 1;
aflow(s2+2,2) = 1;
