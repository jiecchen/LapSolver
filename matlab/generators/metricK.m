function [ a, xy ] = metricK( n )
%METRICK Generates a complete graph weighted by Euclidean distance.
    a = zeros(n,n);
    xy = rand(n,2);
    
    for i = 1:n
       for j = 1:i-1
          a(i,j) = 1 / sqrt( (xy(i,1)-xy(j,1))^2 + (xy(i,2)-xy(j,2))^2 );
          a(j,i) = a(i,j);
       end
    end
    
end

