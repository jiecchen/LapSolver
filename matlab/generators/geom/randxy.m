function xy = randxy(n, type, dimensions)
%RANDXY Returns random points.
    if nargin < 2
        type = 'uniform';
    end
    
    if nargin < 3
        dimensions = 2;
    end
    
    if strcmp(type, 'uniform') 
        xy = rand(n, dimensions);
    elseif strcmp(type, 'normal')
        xy = randn(n, dimensions);
    end
end

