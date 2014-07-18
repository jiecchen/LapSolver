function [ rA ] = incrementalSanitize( rA, k )
%INCREMENTALSANITIZE Summary of this function goes here
%   Detailed explanation goes here

    [u,v,w] = find(rA);
    [~,ind] = sort(w);
    m = length(w);

    for i = 1:m
        x = ind(i);
        rA(u(x), v(x)) = 0;
        rA(v(x), u(x)) = 0;

        if (graphconncomp(rA) >= k)
            break;
        end
        
        w(x);
    end

end

