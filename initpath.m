% script init.m
%
% Created by Dan Spielman on Jun 26, 2007
%
% sets up the correct paths

root = getroot;

paths = strsplit(genpath('matlab'), ':');
paths = paths(1:end-1); % get rid of the terminating ''

for i = 1:length(paths)
    addpath([root, paths(i)]);
end
