function [a, xy] = imgGrid(img, sigma)
%IMGGRID Generates a weighted grid graph from an image.
    dims = size(img);
    
    [a,xy] = grid2(dims(1), dims(2));
    [ai,aj] = find(a);
    
    for e=1:length(ai)
        rgb1 = img(xy(ai(e),1),xy(ai(e),2),:);
        rgb2 = img(xy(aj(e),1),xy(aj(e),2),:);
        delta = rgb1 - rgb2;
        av(e) = exp(-(delta(1)^2 + delta(2)^2 + delta(3)^2) / sigma^2);
    end
    
    a = sparse(ai,aj,av,length(a),length(a));
end

