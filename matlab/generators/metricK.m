function [a, xy] = metricK(n, p)
%METRICK Generates a complete graph with L_p metric weights.
% Default: Euclidean
    a = zeros(n,n);
    xy = rand(n,2);
    
    if nargin < 2
        p = 2;
    end
    
    for i = 1:n
       for j = 1:i-1
          a(i,j) = ( abs(xy(i,1)-xy(j,1))^p + abs(xy(i,2)-xy(j,2))^p )^(1/p);
          a(j,i) = a(i,j);
       end
    end
    
end

