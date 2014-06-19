% A crude comparisson tool between the L matrice outputed by
% Dan's myLDL2 and Serban's LDLDecomposition.

function [a,b] = compLDLmat(graph)
    import lapsolver.algorithms.*;
    import lapsolver.util.*;
    import lapsolver.*;

    n = graph.nv;
    X = rand(1,n);
    
    gvm = GraphVertexRemoval(graph);
    gvmPair = gvm.solve();
    perm = gvmPair.v;
    numRemoved = gvmPair.n;
    g = Graph(GraphUtils.permuteGraph(graph, perm));
    
    la = full(lap(g2a(g)) + diag(X));
    
    ldl = LDLDecomposition(g, X);
    ldlAns = ldl.solve(numRemoved);
    serbanL = full(e2tril(ldlAns.L, n));
    serbanD = full(e2mat(ldlAns.D));
   
    [danL, danD] = myLDL2(la, numRemoved);
    
    keyboard
end
