function gplotComp(a,xy,comp)
% function gplotComp(a,xy,comp)
%
% plot each component of a in a different color.
%
%

minc = min(comp);
maxc = max(comp);

hold on

for i = minc:maxc,
  
  ind = find(comp == i);
  [x,y] = gplot(a(ind,ind),xy(ind,:));
  h =  plot(x,y,'-');
  set(h,'Color',rand(1,3));
  
end

axis off
set(gcf,'Color',[1 1 1]);
