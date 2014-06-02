function [i, j, v] = binLoadIJV(fn)
% function [i, j, v] = binLoadIJV(fn)
%
% exact opposite of binSaveIJV (except weights are floats instead of
% doubles)


fp = fopen(fn,'r');

n = fread(fp, 1, 'int32');
nn = fread(fp, 1, 'int32');

i = fread(fp, nn, 'int32');
i = i + 1;
j = fread(fp, nn, 'int32');
j = j + 1;
v = fread(fp, nn, 'float32');

fclose(fp);

  
