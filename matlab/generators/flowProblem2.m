function aflow = flowProblem2(a,k)
% function aflow = flowProblem2(a,k)
%
% takes as input a graph a (in adj format).
% creates two new nodes, 1 and 2,
% and the problem is to find the maximum flow between them.
%
% each of these nodes are connected to bfs balls containing k
% vertices
%
% variant 1 picks the two nodes to start from at random
% 

n = length(a);

[ai,aj,av] = find(tril(a));

start = ceil(n*rand(1));

pg = yins.Graph(ai,aj);
pg.bfsWalk(start-1);
bfs = pg.bfs + 1;

s1 = bfs(1:k);

last = ceil(n*rand(1));

pg.bfsWalk(last-1);
bfs = pg.bfs + 1;

s2 = bfs(1:k);

aflow((3:n+2),(3:n+2)) = a;
aflow(1,s1+2) = 1;
aflow(s1+2,1) = 1;
aflow(2,s2+2) = 1;
aflow(s2+2,2) = 1;
