function pArray = lowStretchTree(graph)
% function pArray = lowStretchTree(graph)
% 
% computes a low-stretch spanning tree, using 
% Aaron's interface to Sam's code

import graphs.*;

toC = tempname;
fromC = tempname;

binSaveIJV(graph, toC);

host = gethost;
root = getroot;

command = [root, 'barnet/lowstretch2 ', toC, ' ', fromC];

unix(command);

[bi, bj, bv] = binLoadIJV(fromC);

delete(toC);
delete(fromC);

n = length(graph);

b = sparse(bi,bj,bv);
b(n,n) = 0;
b = b + b';
[bi,bj,bv] = find(tril(b));
pg = PreGraph(bi,bj,bv);
pArray = pg.treeToArray;

