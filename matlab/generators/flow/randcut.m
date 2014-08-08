function ac = randcut(a,s,t)
%RANDCUT Finds the adjacency matrix of a random cut in a graph.
    import lapsolver.algorithms.*;

    n = length(a);
    g = a2g(a);
    
    if nargin == 3
        cut = RandomCut.compute(g,s,t);
    else
        cut = RandomCut.compute(g);
    end

    ac = sparse(double(cut.u+1), double(cut.v+1), cut.weight, n, n);
    ac = ac + ac';
end

