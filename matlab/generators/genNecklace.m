function ab = genNecklace(a, b)
% function ab = genNecklace(a, b)
%
% this is the "necklace" product of a and b
% it makes numverts(b) copies (randpermed) of a,
%
% and then connects them according to b,
% with random hookups.
%

na = length(a);
nb = length(b);

aneck = sparse(na*nb,na*nb);

for i = 1:nb,
    p = randperm(na);
    s = (1:na)+(na*(i-1));
    ab(s,s) = a(p,p);
end

[bi,bj] = find(tril(b));

E = zeros(length(bi),2);

for i = 1:length(bi),
    x1 = na*(bi(i)-1) + ceil(na*rand(1));
    x2 = na*(bj(i)-1) + ceil(na*rand(1));
    E(i,:) = [x1,x2];
end

ab2 = sparse(E(:,1),E(:,2),1,na*nb,na*nb);
ab = ab + ab2 + ab2';

