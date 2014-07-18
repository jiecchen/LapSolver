function nA = mtNormalize(A, rA, f, e, k, eps)
%MTNORMALIZE Performs the Miller-Tolliver normalization step.
%   Detailed explanation goes here
    nA = rA;
    alpha = 1;
    threshold = mtNormFA(A, A, f, e, k) + eps;
    while mtNormFA(nA, A, f, e, k) > threshold
        alpha = alpha / 2;
        nA = (1-alpha)*A + alpha*nA;
    end
end

