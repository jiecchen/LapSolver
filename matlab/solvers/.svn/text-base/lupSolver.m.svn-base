function x = lupSolver(l,u,p,b)
% function x = lupSolver(l,u,p,b)
% function f = lupSolver(l,u,p)
%
% given an lup factorization, applies it to solve for x
% if b is not given, return a function that does this.
%
% for example, 
% [l,u,p] = lu(la);
% x = lupSolver(l,u,p,b);
%


default('b',[]);

% if p is a perm matrix, turn it into a vector
if (prod(size(p))>length(p)),
    [p,jnk] = find(p);
end

if isempty(b)
    f = @(b)(internal(l,u,p,b));
    x = f;
else
    x = internal(l,u,p,b);
end

end % main function

function x = internal(l,u,p,b)

  bp = b(p);
  yp = l \ bp;
  xp = u \ yp;
  x(p) = xp;
  x = x(:);

end
  