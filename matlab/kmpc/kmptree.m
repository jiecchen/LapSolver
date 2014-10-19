function [ at ] = kmptree( a )
%KMPTREE
nameIn = strcat(tempname, '_in.adj');
nameOut = strcat(tempname, '_out.adj');
ijvwrite(nameIn, a);
system(sprintf('./kmp-c/build/lsst < %s > %s', nameIn, nameOut));
at = ijvread(nameOut);
end
