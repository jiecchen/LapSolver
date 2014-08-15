function a = ijvread(filename)
%IJVREAD Reads a graph in C adjacency list form.
    ve = dlmread(filename, ' ', [0 0 0 1]);
    nv = ve(1);
    ne = ve(2);
    
    ijv = dlmread(filename, ' ', 1, 0);
    
    if max(max(ijv(:,1:2))) ~= nv-1
        fprintf('warning: inconsistent vertex count');
    end
    
    if size(ijv,1) ~= ne
        fprintf('warning: inconsistent edge count');
    end
    
    a = sparse(ijv(:,1)+1, ijv(:,2)+1, ijv(:,3), nv, nv);
    a = a + a';
    
    if graphconncomp(a) ~= 1
        fprintf('warning: graph not connected');
    end
end

