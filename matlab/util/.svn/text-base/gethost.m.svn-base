function host = gethost()
% function host = gethost()

host = getenv('HOST');
if isempty(host)
  host = getenv('HOSTNAME');
end
if isempty(host)
  host = 'HOST';
end

