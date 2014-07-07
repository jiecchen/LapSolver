function x = maxflow(A, s, t, c) 
	% function x = maxflow(A, s, t, c)

	[A, c, n, m, u_t] = ProcessGraph(A, s, t, c);

	U = max(c(1 : m - 1));

	ep = 0.5;
	ep_flow = ep^2 / (64 * m^2 * n^2 * U^3);

	A_new = [A, sparse(n, m), sparse(n, m), speye(n), -speye(n);
		     speye(m), speye(m), -speye(m), sparse(m, n), sparse(m, n)];

	b = full([zeros(n, 1); c]);

	k = 4 * U / ep_flow;
	height = m + 2 * n;
	c = [-u_t; zeros(m, 1); k * ones(height, 1)];

	lmin = 2;
	T = (n * U + 1) * 4 * U / ep_flow + 1;

	k = -2 * U / ep_flow;
	y = [zeros(n, 1); k * ones(m, 1)];

	%%%%%%%%%%%%%%% interior point %%%%%%%%%%%%%%%%%%%%

	M = 3 * m + 2 * n;

	s = c - [k * ones(2 * m, 1); -k * ones(m, 1); zeros(2 * n, 1)];
	s_gap = M;

	s_inv = s.^(-1);
	b_ = full(A_new) * s_inv;
	z = b_' * y - s_gap;

	bound = 40 * T * M * norm(b_) / sqrt(lmin);

	iter1 = 0;

	while s_gap < bound
		z = z - s_gap / log(m);
		% z = z - s_gap / (10 * sqrt(M));

		s1 = c - [A' * y(1 : n) + y(n + 1 : end); y(n + 1 : end); -y(n + 1 : end); y(1 : n); -y(1 : n)];
		s2 = -z + b_' * y;
		s2 = s2 * ones(M, 1);
		s = [s1; s2];

		s_inv = s.^(-1);
		s_inv2 = s_inv.^(2);

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		SA = spdiags(s_inv2(1 : m), 0, m, m);
		SAAt = SA * A';
		W = [[A * SAAt + spdiags(s_inv2(3 * m + 1 : 3 * m + n) + s_inv2(3 * m + n + 1 : M), 0, n, n), A * SA];
			 [SAAt, spdiags(s_inv2(1 : m) + s_inv2(m + 1 : 2 * m) + s_inv2(2 * m + 1 : 3 * m), 0, m, m)]];

		r = [-A * s_inv(1 : m) - s_inv(3 * m + 1 : 3 * m + n) + s_inv(3 * m + n + 1 : M);
			 -s_inv(1 : m) - s_inv(m + 1 : 2 * m) + s_inv(2 * m + 1 : 3 * m)];
		r = r + b_ * sum(s_inv(M + 1 : end)); %combine sum term

		sigma = sum(s_inv2(M + 1 : end));

		% [result1, garbage] = pcg(W, r, 1 / (20 * (sqrt(M) + 1)), 40);
		% [result2, garbage] = pcg(W, b_, 1 / (20 * (sqrt(M) + 1)), 40);
		[result1, garbage] = cmgSolver(W, r);
		[result2, garbage] = cmgSolver(W, b_);


		% result1 = pcg(W, r, 1 / (20 * (sqrt(M) + 1)), 40);
		% result2 = pcg(W, b_, 1 / (20 * (sqrt(M) + 1)), 40);

		result = result1 - (sigma * (result2 * (b_' * result1)) / (1 + sigma * (b_' * result2)));

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


		y = y + result * (1 - 1 / (20 * (sqrt(M) + 1)));
		s_gap = b_' * y - z;

		iter1 = iter1 + 1;
	end

	iterations_of_unshift = iter1;

	z = b' * y - 40 * T * M * norm(b) / sqrt(lmin);
	s_gap = b' * y - z;

	iter2 = 0;

	while s_gap > ep / 3
		z = z + s_gap / log(m); 
		% z = z + s_gap / (10 * sqrt(M));

		s1 = c - [A' * y(1 : n) + y(n + 1 : end); y(n + 1 : end); -y(n + 1 : end); y(1 : n); -y(1 : n)];
		s2 = -z + b' * y;
		s2 = s2 * ones(M, 1);
		s = [s1; s2];
		s_inv = s.^(-1);
		s_inv2 = s_inv.^(2);

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		SA = spdiags(s_inv2(1 : m), 0, m, m);
		SAAt = SA * A';
		W = [[A * SAAt + spdiags(s_inv2(3 * m + 1 : 3 * m + n) + s_inv2(3 * m + n + 1 : M), 0, n, n), A * SA];
			 [SAAt, spdiags(s_inv2(1 : m) + s_inv2(m + 1 : 2 * m) + s_inv2(2 * m + 1 : 3 * m), 0, m, m)]];

		r = [-A * s_inv(1 : m) - s_inv(3 * m + 1 : 3 * m + n) + s_inv(3 * m + n + 1 : M);
			 -s_inv(1 : m) - s_inv(m + 1 : 2 * m) + s_inv(2 * m + 1 : 3 * m)];
		r = r + b * sum(s_inv(M + 1 : end)); %combine sum term

		sigma = sum(s_inv2(M + 1 : end));


		% [result1, garbage] = pcg(W, r, 1 / (20 * (sqrt(M) + 1)), 40);
		% [result2, garbage] = pcg(W, b, 1 / (20 * (sqrt(M) + 1)), 40);
		[result1, garbage] = cmgSolver(W, r);
		[result2, garbage] = cmgSolver(W, b);
		% result1 = pcg(W, r, 1 / (20 * (sqrt(M) + 1)), 40);
		% result2 = pcg(W, b, 1 / (20 * (sqrt(M) + 1)), 40);

		result = result1 - (sigma * (result2 * (b' * result1)) / (1 + sigma * (b' * result2)));

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		y = y + result * (1 - 1 / (20 * (sqrt(M) + 1)));
		s_gap = b' * y - z;

		iter2 = iter2 + 1;
	end

	iterations_of_shift = iter2;

	% s = c - [A' * y(1 : n) + y(n + 1 : end); y(n + 1 : end); -y(n + 1 : end); y(1 : n); -y(1 : n)];
	% s_inv = s.^(-1);
	% s_inv2 = s_inv.^(2);

	s1 = c - [A' * y(1 : n) + y(n + 1 : end); y(n + 1 : end); -y(n + 1 : end); y(1 : n); -y(1 : n)];
	s2 = b' * y - z;
	s2 = s2 * ones(M, 1);
	s = [s1; s2];
	s_inv = s.^(-1);
	s_inv2 = s_inv.^(2);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	SA = spdiags(s_inv2(1 : m), 0, m, m);
	SAAt = SA * A';
	W = [[A * SAAt + spdiags(s_inv2(3 * m + 1 : 3 * m + n) + s_inv2(3 * m + n + 1 : M), 0, n, n), A * SA];
		 [SAAt, spdiags(s_inv2(1 : m) + s_inv2(m + 1 : 2 * m) + s_inv2(2 * m + 1 : 3 * m), 0, m, m)]];

	r = [A * s_inv(1 : m) + s_inv(3 * m + 1 : 3 * m + n) - s_inv(3 * m + n + 1 : M);
		 s_inv(1 : m) + s_inv(m + 1 : 2 * m) - s_inv(2 * m + 1 : 3 * m)];
	r = r - b * sum(s_inv(M + 1 : end));

	% r = [-A * s_inv(1 : m) - s_inv(3 * m + 1 : 3 * m + n) + s_inv(3 * m + n + 1 : M);
	% 	 -s_inv(1 : m) - s_inv(m + 1 : 2 * m) + s_inv(2 * m + 1 : 3 * m)];
	% r = r + b * sum(s_inv(M + 1 : end)); %combine sum term

	sigma = sum(s_inv2(M + 1 : end));


	% [result1, garbage] = pcg(W, r, .1, 40);
	% [result2, garbage] = pcg(W, b, .1, 40);
	[result1, garbage] = cmgSolver(W, r);
	[result2, garbage] = cmgSolver(W, b);
	% result1 = pcg(W, r, .1, 40);
	% result2 = pcg(W, b, .1, 40);

	result = result1 - (sigma * (result2 * (b' * result1)) / (1 + sigma * (b' * result2)));

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	x = s_inv(1 : M) - spdiags(s_inv2(1 : M), 0, M, M) * A_new' * result; %fix this later
	x = x / M;
	x = x / ((1 / s_gap) + (b' * result / (s_gap^(2))));

	x = round(norm(x(m)));
end

function [A, c, n, m, u_t] = ProcessGraph(A, s, t, c)
	% makes sure there is exactly one edge leading into t,
	% changes A into an incidence matrix from a graph object
	% modifies the capacity vector to accomodate a new edge
	% removes the rows corresponding to s and t and adjusts n and m


	A = a2u(A);
	n = size(A, 1);
	m = size(A, 2);

	% full(A)

	new_cost = A(t, :) * c;

	A = [A, zeros(n, 1)]; % add new column for edge from old t to new
	m = m + 1;
	A(t, m) = -1;

	A(s, :) = [];
	A = sparse(A);
	n = n - 1;

	c = [c; new_cost]; % final edge into t cannot be min cut

	u_t = zeros(m, 1); % zero vector with one corresponding to edge into t
	u_t(m, 1) = 1;
	u_t = sparse(u_t);
end




