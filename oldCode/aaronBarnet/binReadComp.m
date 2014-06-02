function comp = binReadComp(fn)
% function binReadComp(fn)
% read a binary file listing components.
% first int gives number of nodes
% rest of ints give comps for corresponding verts


fp = fopen(fn,'r');

n = fread(fp, 1, 'int32');

comp = fread(fp, n, 'int32');


fclose(fp);

  
