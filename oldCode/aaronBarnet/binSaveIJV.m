function binSaveIJV(a,fn)
% function binSaveIJV(a,fn)
%
% save graphs "a" as file fn in a binary format,
% to be read by binReadIJV
%
% first comes n, as an int32 (signed),
% then nzz,  as an int32 (signed),
% then i, as a vector, in int32 (signed), based a 1
% then j, as a vector, in int32 (signed), based a 1
% then v, as a vector, in 64 bit floats


fp = fopen(fn,'w');

[i,j,v] = find(tril(a));

fwrite(fp, length(a), 'int32');
fwrite(fp, length(i), 'int32');

fwrite(fp, i, 'int32');
fwrite(fp, j, 'int32');
fwrite(fp, v, 'float64');


fclose(fp);

  
