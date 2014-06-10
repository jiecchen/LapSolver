function gplot3(a, xyz)
%GPLOT3 Plots a graph in 3D, given adjacency matrix and 3-by-n coords
%matrix
    [i, j] = find(triu(a));
    X = permute(cat(3, xyz(i,:), xyz(j,:)), [3 1 2]);
    plot3(X(:,:,1),X(:,:,2),X(:,:,3));
end

