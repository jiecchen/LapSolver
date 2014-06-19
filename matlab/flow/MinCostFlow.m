function x = MinCostFlow(A, s, t, c, q, n, m, ep)
	% function x = MinCostFlow(A, s, t, c, q, m, n, ep)
	% A is an nxm incidence matrix
	% s and t are the start and end of the min cost flow t > s
		% there must be exactly one edge into t
	% c is a length m vector of edge capacities
	% q is a length m vector of edge costs
	% m is the number of edges
	% n is the number of vertices
	% ep is the additive error tolerance

	maxflow = MaxFlow(A, s, t, c, n, m, ep);
	F = A * maxflow;
	F = F(t, 1);

	A(s, :) = []; % careful with this...
	t = t - 1;
	et = zeros(n-1, 1);
	et(t, 1) = 1;

	U = max([max(c), max(q)]);
	%epflow = ep^2 / (64 * m^2 * n^2 * U^3);
	epflow = ep;

	Aip = [A, zeros(n-1, m), zeros(n-1, m), eye(n-1, n-1), -1 * eye(n-1, n-1);
		   eye(m, m), eye(m, m), -1 * eye(m, m), zeros(m, n-1), zeros(m, n-1)];
	bip = [F * et; c];

	p = 4 * m * U^2 / epflow;
	cip = [q; zeros(m, 1); p * ones(m, 1); p * ones(n-1, 1); p * ones(n-1, 1)];
	
	lmin = 2;
	T = 5; % ???

	p = -m * U^2 / epflow;
	yip = [zeros(n-1, 1); p * ones(m, 1)];
	
	xip = InteriorPoint(Aip, bip, cip, lmin, T, yip, ep);

	x = xip(1:m);
end






