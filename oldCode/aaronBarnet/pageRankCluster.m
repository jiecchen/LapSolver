function comp = pageRankCluster(a, clusSize, alpha, epsilon_divider)
% function comp = pageRankCluster(a, clusSize, alpha, epsilon_divider)
%
% run Aaron's graphcluster code
%
% for example, try this:
%
% [a,bdry,xy] = grid2(19);
% comp = pageRankCluster(a,12, .125, 4.0);
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

%[i, j, v] = find(tril(a));

%saveijv(i, j, v, length(a), graphfile);
binSaveIJV(a, graphfile);
root = getroot;


cmd = [root, 'barnet/graphclusterpr ', num2str(clusSize), ' ', num2str(alpha), ' ', num2str(epsilon_divider),' ', graphfile, ...
       ' ' compfile];

system(cmd);

%comp = compfile
comp = binReadComp(compfile);

delete(graphfile);
delete(compfile);
