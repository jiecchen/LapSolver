A = full(grid2(10));
n = size(A, 1);
D = diag(sum(A));
L = D - A;
I = diag(ones(n,1));
nL = (I+D)^(-1/2)*(I+L)*(I+D)^(-1/2);

numItr = 10;
B = cell(numItr, 1);
R = I - nL;
for i=1:numItr
    B{i} = I + R./2;
    R = I - (I + R./2)*(I - R)*(I + R./2);
end

res = I;
for i=1:numItr
    res = res * B{i};
end

norm(nL^(-0.5)-res)