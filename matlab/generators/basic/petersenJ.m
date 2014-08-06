function [ a, coords ] = petersenJ( n, k )
% function [a,corrds] = petersenJ( n, k )
%
% a G(n,k) Generalized Petersen Graph
%
    import lapsolver.generators.PetersenGraph;
   
    gen = PetersenGraph(n,k);
    a = g2a(gen.generateGraph);
    
    ns = 2*pi*(0:(n-1))/n;
    circle = [cos(ns); sin(ns)];
    
    coords = [circle, circle/2]';
end

