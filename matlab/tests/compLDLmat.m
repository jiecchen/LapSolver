% A crude comparisson tool between the L matrice outputed by
% Dan's myLDL2 and Serban's LDLDecomposition.

function compLDLmat(graph)
    import lapsolver.algorithms.*;
    import lapsolver.util.*;
    import lapsolver.*;

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
    serbanL = full(e2tril(ldlAns.L, n));
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
        keyboard
    else
        fprintf('%d %d GG\n', graph.nv, graph.ne);
    end
    
    keyboard
end
