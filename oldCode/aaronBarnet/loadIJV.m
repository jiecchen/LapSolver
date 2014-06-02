function a = loadIJV(fn)
% function a = loadIJV(fn)
%
% load graphs "a" from file fn in the following format:
%  first line is 
% n m  -- where n is number of vertices and m is number of edges
%
% this is followed by m lines, each of the form
% from to weight
%
% returns symmetrized version


fp = fopen(fn,'r');


nv = fscanf(fp,'%d',1)
ne = fscanf(fp,'%d',1)

ai = zeros(ne,1);
aj = zeros(ne,1);
av = zeros(ne,1);


for x = 1:ne,
  [trip] = fscanf(fp,'%d %d %f\n',3);
  ai(x) = trip(1);
  aj(x) = trip(2);
  av(x) = trip(3);
end

fclose(fp);

a = sparse(ai,aj,av,nv,nv);
a = a + a';

  
