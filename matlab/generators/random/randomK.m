function a = randomK( n )
%RANDOMK Generates a clique with random weights.
    a = tril(rand(n,n),-1);
    a = a + a';
end
