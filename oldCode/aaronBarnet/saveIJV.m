function saveIJV(a,fn)
% function saveIJV(a,fn)
%
% save graphs "a" as file fn in the following format:
%  first line is 
% n m  -- where n is number of vertices and m is number of edges
%
% this is followed by m lines, each of the form
% from to weight


fp = fopen(fn,'w');
[ind,jnk] = find(tril(a));

fprintf(fp,'%d %d\n',length(a),length(ind));

[i,j,v] = find(tril(a));
for x = 1:length(i),
  fprintf(fp,'%d %d %g\n',i(x),j(x),v(x));
end

fclose(fp);

  
