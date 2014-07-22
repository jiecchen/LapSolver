function [ partA ] = smartCluster(A, k)
%SMARTCLUSTER Uses a maximum spanning tree to form k clusters for a given
%             graph

    n = length(A);
    [u, v, w] = find(A);
    
    % As all the edge costs in A are less than 1, I will use minimum 
    % spanning tree on 1/edge costs
    
    m = length(w);
    for i = 1:m
        w(i) = -w(i) + 10;
    end
    
    A = sparse(u, v, w);
    
    [t, p] = graphminspantree(A);
    
    % I will remove the largest k-1 edges from the tree
    
    [ut, vt, wt] = find(t);
    [~, ind] = sort(wt);
    
    for i = (n-k+1):(n-1)
        wt(ind(i)) = 0;
    end
    
    for i = 1:(n-1)
        if ut(i) > vt(i)
            aux = ut(i);
            ut(i) = vt(i);
            vt(i) = aux;
        end
    end
    
    t = sparse(ut, vt, wt, n, n);
    t = t + t';
    
    % I will create the answer, broken in k components
    
    [s, c] = graphconncomp(t);
    
    ufin = u;
    vfin = v;
    wfin = w;
    
    for i = 1:m
        if c(u(i)) ~= c(v(i))
            wfin(i) = 0;
        end
    end
    
    partA = sparse(ufin, vfin, wfin);
    
    %keyboard;
end

