function ar = perturb(a, eps)
%PERTURB Perturbs weights by a small additive term.
    if nargin < 2
        eps = 1e-6;
    end

    ar = sprand(tril(a,-1)) * eps;
    ar = a + ar + ar';
end

