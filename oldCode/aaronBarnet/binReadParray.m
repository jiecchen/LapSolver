function pa = binReadParray(fn)
% function binReadComp(fn)
% read a binary file listing components.
% first int gives number of nodes
% rest of ints give  parent pointers for corresp vertices
% need to add 1 to index for matlab


fp = fopen(fn,'r');

n = fread(fp, 1, 'int32');

pa = fread(fp, n, 'int32');
pa = pa + 1;

fclose(fp);

  
