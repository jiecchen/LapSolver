% A crude comparisson tool between the L matricx outputed by
% Dan's myLDL2 and Serban's LDLDecomposition.

function compLDLmat(graph)
    import lapsolver.algorithms.*;
    import lapsolver.util.*;
    import lapsolver.*;

    % multiply edge waits by some random values
    mat = g2a(graph);
    [u,v,w] = find(mat);
    for i = 1:length(w)
        w(i) = w(i) * rand(1,1);
    end
    mat = sparse(u,v,w);
    graph = a2g(mat);
    
    n = graph.nv;
%    lagraph = full(lap(g2a(graph)));
    X = rand(1,n);
    
    gvm = GraphVertexRemoval(graph);
    gvmPair = gvm.solve();
    perm = gvmPair.permutation;
    numRemoved = gvmPair.numRemoved;
    g = Graph(GraphUtils.permuteGraph(graph, perm));
    
    ldl = LDLDecomposition(g, X);
    ldlAns = ldl.solve(numRemoved);
    serbanL = full(e2mat(ldlAns.L.L));
    serbanD = full(e2mat(ldlAns.D));
    
    la = full(lap(g2a(g)) + diag(X));
   
    [danL, danD] = myLDL2(la, numRemoved);
    
    eps = 10^(-10);
    valueL = max(max(abs(serbanL - danL)));
    valueD = max(max(abs(serbanD - danD)));
    
    valueL = 0;
    valueD = 0;
    if valueL > eps || valueD > eps
        fprintf('%d %d WA\n', graph.nv, graph.ne);
    else
        fprintf('%d %d GG\n', graph.nv, graph.ne);
    end
end
