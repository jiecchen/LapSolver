function x = MaxFlow(A, s, t, c, m, n, ep)
	% function x = MaxFlow(A, s, t, c, m, n, ep)
	% a nxm directed graph A
	% indices s and t for a max flow from s to t
	% NOTE: there must be exactly one edge directed into t
	% length m vector c corresponding to positive integer edge capacities
	% m := number of edges
	% n := number of vertices

	ut = A(t, :);
	ut = ut';
	A(s, :) = [];
	A(t-1, :) = [];

	U = max(c);
	epflow = ep^2 / (64 * m^2 * n^2 * U^3);

	Aip = [A, zeros(n-2, m), zeros(n-2, m), eye(n-2, n-2), -1 * eye(n-2, n-2);
		   eye(m, m), eye(m, m), -1 * eye(m, m), zeros(m, n-2), zeros(m, n-2)];
	bip = [zeros(n-2, 1); c];
	p = 4 * U / epflow;
	cip = [ut; zeros(m, 1); p * ones(m, 1); p * ones(n-2, 1); p * ones(n-2, 1)];
	lmin = 2;
	T = 1000; % ???
	p = -2 * U / epflow;
	yip = [zeros(n-2, 1); p * ones(m, 1)];

	xip = InteriorPoint(Aip, bip, cip, lmin, T, yip, ep);
	x = xip(1:m);
end
