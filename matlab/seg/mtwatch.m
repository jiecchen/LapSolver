function [ar, w] = mtwatch(a, k, iters)
%MTWATCH Watch the evolution of the edge weights in Miller-Tolliver.
    w = zeros(iters, nnz(a));
    ar = a;
    for i = 1:iters
        [~,~,av] = find(ar);
        w(i,:) = av;
        ar = mtGraph(ar, k, 1e-3);
    end
end

