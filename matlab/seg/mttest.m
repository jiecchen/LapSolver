function [ar, colors, a] = mttest(path, k, tol, edgetol, sigma)
%MTTEST Test for Miller-Tolliver image segmentation.
    image = double(imread(path));
    dims = size(image);
    
    if nargin < 5
        sigma = 70;
    end
    
    a = imgGrid(image, sigma);
    ar = mtGraph(a, k, tol, edgetol);
    
    [~, colors] = graphconncomp(ar);
    colors = reshape(colors, dims(1), dims(2));
end

