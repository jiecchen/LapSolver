function root = getroot()
% function root = getroot()
%
% get linsolve root directory
%
% if LINSOLVE is one of your environment variables, will make this
% the root
% for example, in .bashrc have:
%  export LINSOLVE=/export/home/spielman/rep/linsolve/
%
% if LINSOLVE is not an env variable, will use the current directory

root = getenv('LINSOLVE');
if (isempty(root)),
  root = pwd;
end

if (isempty(root))
  root = pwd;
  display(['making the current directory the root']);
end

if (root(end) ~= '/')
  root = [root, '/'];
end
