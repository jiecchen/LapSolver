function a = eyeFlowGraph(k)
% function a = eyeFlowGraph(k)
%
% a graph that is bad for some flow algorithms:
% should take k iterations for interior point or linsolve based
% algorithms
%
% will have k^3 nodes
%
% the problem is to flow from 1 to 2
%
% still off by a little

if (k <= 1)
    e = [1 3; 2 3; 1 4; 2 4];
    a = sparse(e(:,1),e(:,2),1,4,4);
    a = a + a';
    return
end

a0 = eyeFlowGraph(k-1);
n0 = length(a0);

% bundle = diag(sparse(ones(1,k*(k-1))),k);
if (k > 2),
  bundle = diag(sparse(ones(1,k*(k-2))),k);
  bundle = bundle+bundle';
else
    bundle = sparse(2,2);
end


nb = length(bundle);

a1 = [a0, zeros(n0,2*nb); 
      zeros(nb,n0), bundle, zeros(nb,nb);
      zeros(nb,n0+nb), bundle];
n1 = length(a1);

% so, recursive part has terminals 3 and 4
a = [zeros(2,n1+2); zeros(n1,2), a1];

sb1 = 2+n0;
sb2 = 2+n0+nb;

a(1,4) = 1;
a(4,1) = 1;

% attach bundle1 to nodes 1 and 3
a(1,sb1+[1:k]) = 1;
a(sb1+[1:k],1) = 1;
a(3,sb1+nb-k+[1:k]) = 1;
a(sb1+nb-k+[1:k],3) = 1;


a(2,3) = 1;
a(3,2) = 1;

% attach bundle2 to nodes 2 and 4
a([4],sb2+[1:k]) = 1;
a(sb2+[1:k],[4]) = 1;
a([2],sb2+nb-k+[1:k]) = 1;
a(sb2+nb-k+[1:k],[2]) = 1;



