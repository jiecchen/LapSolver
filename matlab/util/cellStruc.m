function y = cellStruc(x)
% function y = cellStruc(x)
%
% if x is a cell array of strucs,
% y is a struc of cell arrays of the same size

l = length(x);

fn = fieldnames(x{1});

y = [];

for f = 1:length(fn),
  
  dat = zeros(1,l);
  for i = 1:l,
    gf = getfield(x{i},fn{f});
    if (isstr(gf))
      dat(i) = NaN;
    else
      dat(i) = gf;
    end
    
  end
  
  y = setfield(y,fn{f},dat);
end
