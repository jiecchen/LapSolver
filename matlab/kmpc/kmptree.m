function [ at ] = kmptree( a )
%KMPTREE
nameIn = strcat(tempname, '_in.adj');
nameOut = strcat(tempname, '_out.adj');
ijvwrite(nameIn, a);
system(sprintf('source /opt/intel/composerxe/bin/compilervars.sh intel64 && ./kmp-c/build/lsst < %s > %s', nameIn, nameOut));
at = ijvread(nameOut);
end
