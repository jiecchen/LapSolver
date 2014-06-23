function x = MaxFlow(A, s, t, c)
	% function x = MaxFlow(A, s, t, c, ep)
	% A is a graph object
	% s and t are the beginning and end of the flow
	% c is a length m vector of edge capacities

	[A, c, n, m, u_t] = ProcessGraph(A, s, t, c);

	U = max(c(1:m-1));

	ep = 0.5;
	ep_flow = ep^2 / (64 * m^2 * n^2 * U^3);

	A_new = [A, zeros(n, m), zeros(n, m), eye(n), -eye(n);
		     eye(m), eye(m), -eye(m), zeros(m, n), zeros(m, n)];

	b_new = [zeros(n, 1); c];

	k = 4 * U / ep_flow;
	height = m + 2 * n;
	c_new = [-u_t; zeros(m, 1); k * ones(height, 1)];

	lmin = 2;
	T = (n * U + 1) * 4 * U / ep_flow + 1;

	k = -2 * U / ep_flow;
	y_new = [zeros(n, 1); k * ones(m, 1)];

	% x = InteriorPoint(A_new, b_new, c_new, lmin, T, y_new, ep);
	x = ip(A_new, b_new, c_new, y_new, T, lmin, ep);

	x = x(1:m);
	x = u_t' * x;
end

function [A, c, n, m, u_t] = ProcessGraph(A, s, t, c)
	% makes sure there is exactly one edge leading into t,
	% changes A into an incidence matrix from a graph object
	% modifies the capacity vector to accomodate a new edge
	% removes the rows corresponding to s and t and adjusts n and m


	A = a2u(A);
	n = size(A, 1);
	m = size(A, 2);

	A = [A, zeros(n, 1)]; % add new column for edge from old t to new
	m = m + 1;
	A(t, m) = -1;

	A(s, :) = [];
	n = n - 1;

	c = [c; m * max(c)]; % final edge into t cannot be min cut

	u_t = zeros(m, 1); % zero vector with one corresponding to edge into t
	u_t(m, 1) = 1;
end



















