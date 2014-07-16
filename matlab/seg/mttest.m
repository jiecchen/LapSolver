function [ar, colors, a] = mttest(path, k, tol, edgetol)
%MTTEST Summary of this function goes here
    image = double(imread(path));
    dims = size(image);
    
    a = imgGrid(image, 70);
    ar = mtGraph(a, k, tol, edgetol);
    
    [~, colors] = graphconncomp(ar);
    colors = reshape(colors, dims(1), dims(2));
end

