function comp = greedyCluster(a, clusSize)
% function comp = greedyCluster(a, clusSize)
%
% run Aaron's graphcluster code
%
% for example, try this:
%
% [a,bdry,xy] = grid2(19);
% comp = greedyCluster(a,12);
% gplotComp(a,xy,comp)
%
% and compare to Dan's old java code by
%
% import graphs.*;
% import cluster.*;
% [ai,aj,av] = find(tril(a));
% pg = PreGraph(ai,aj,av);
% cl = ClusterV(pg,2);
% cl.cluster2(0, 12);
% comp = cl.comp;
% gplotComp(a,xy,comp)
%


graphfile = tempname;
compfile = tempname;

saveIJV(a,graphfile);

root = getroot;

cmd = [root, 'barnet/graphcluster ', num2str(clusSize), ' ', graphfile, ...
       ' ' compfile];

system(cmd);

comp = load(compfile);

delete(graphfile);
delete(compfile);
