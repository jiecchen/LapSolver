function x = copyStruct(x,y)
% function x = copyStruct(x,y)
%
% copy all of the fields of y onto x

nm = fieldnames(y);
for i = 1:length(nm),
  x = setfield(x,nm{i},getfield(y,nm{i}));
end
