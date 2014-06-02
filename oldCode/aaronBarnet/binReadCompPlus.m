function [comp, used, rej] = binReadComp(fn)
% function binReadComp(fn)
% read a binary file listing components.
% first int gives number of nodes
% rest of ints give comps for corresponding verts

fp = fopen(fn,'r');

n = fread(fp, 1, 'int32');
comp = fread(fp, n, 'int32');

n = fread(fp, 1, 'int32');
used_size = fread(fp, n, 'int32');
%n = fread(fp, 1, 'int32');
used_conduct = fread(fp, n, 'float64');
used = [used_size, used_conduct];


n = fread(fp, 1, 'int32');
rej_size = fread(fp, n, 'int32');
%n = fread(fp, 1, 'int32');
rej_conduct = fread(fp, n, 'float64');
rej = [rej_size, rej_conduct];

fclose(fp);