function x = maxflow(A, s, t, c)
	% function x = maxflow1(A, s, t, c)
	[A, c, n, m, u_t] = ProcessGraph(A, s, t, c);

	U = max(c(1 : m - 1));

	ep = 0.5;
	ep1 = 1e-6;
	ep_flow = ep^2 / (64 * m^2 * n^2 * U^3);

	b = [zeros(n, 1); c];

	k = 4 * U / ep_flow;
	height = m + 2 * n;
	c = [-u_t; zeros(m, 1); k * ones(height, 1)];

	lmin = 2;
	T = (n * U + 1) * 4 * U / ep_flow + 1;

	k = -2 * U / ep_flow;
	y = [zeros(n, 1); k * ones(m, 1)];

	%%%%%%%%%%%%%%% interior point %%%%%%%%%%%%%%%%%%%%

	M = 3 * m + 2 * n;

	%%% FindCentralPath() %%%
	s = c - [k * ones(2 * m, 1); -k * ones(m, 1); zeros(2 * n, 1)];
	s_gap = M;

	s_inv = s.^(-1);
	b_ = [A * s_inv(1 : m) + s_inv(3 * m + 1 : 3 * m + n) - s_inv(3 * m + n + 1 : M);
		  s_inv(1 : m) + s_inv(m + 1 : 2 * m) - s_inv(2 * m + 1 : 3 * m)];
	z = b_' * y - s_gap;

	bound = 40 * T * M * norm(b_) / sqrt(lmin);
	iter1 = 0;

	while s_gap < bound
		%%% Unshift() %%%
		z = z - s_gap / log(m);
		% z = z - s_gap / (10 * sqrt(M));
			%%% NewtonStep() %%%
			s = c - [A' * y(1 : n) + y(n + 1 : end); y(n + 1 : end); -y(n + 1 : end); y(1 : n); -y(1 : n)];
				%%% Solve(A*S_inv2*A', -As_inv) %%%
				s_inv = s.^(-1);
				r = [-A * s_inv(1 : m) - s_inv(3 * m + 1 : 3 * m + n) + s_inv(3 * m + n + 1 : end);
					 -s_inv(1 : m) - s_inv(m + 1 : 2 * m) + s_inv(2 * m + 1 : 3 * m)];
				r = r + (M / (b_' * y - z)) * b_;

				s_inv2 = s.^(-2);
				D1 = spdiags(s_inv2(1 : m), 0, m, m);
				D2 = spdiags(s_inv2(3 * m + 1 : 3 * m + n) + s_inv2(3 * m + n + 1 : M), 0, n, n);
				d13 = s_inv2(1 : m) + s_inv2(m + 1 : 2 * m) + s_inv2(2 * m + 1 : 3 * m);
				d13 = spdiags(d13.^(-1), 0, m, m);

				Ms = A * D1 * A' + D2 - A * D1 * d13 * D1 * A';
				Q = A * D1 * d13;
				v = sqrt(M) / (b_' * y - z) * b_;

				% x1 = cmgSolver(Ms, r(1 : n) - Q * r(n + 1 : end));
				[x1, garbage] = pcg(Ms, r(1 : n) - Q * r(n + 1 : end), ep1);
				x1 = [x1; d13 * (r(n + 1 : end) - D1 * A' * x1)];

				% x2 = cmgSolver(Ms, v(1 : n) - Q * v(n + 1 : end));
				[x2, garbage] = pcg(Ms, v(1 : n) - Q * v(n + 1 : end), ep1);
				x2 = [x2; d13 * (v(n + 1 : end) - D1 * A' * x2)];

				x = x1 - (x2 * (x2' * r)) / (1 + v' * x2);
				%%% End Solve() %%%
			y = y + (1 - 1 / (20 * (sqrt(M) + 1))) * x;
			%%% End NewtonStep() %%%
		s_gap = b_' * y - z;	
		%%% End Unshift() %%%

		iter1 = iter1 + 1;
	end
	iterations_of_unshift = iter1

	z = b' * y - 40 * T * M * norm(b) / sqrt(lmin);
	s_gap = b' * y - z;
	%%% End FindCentralPath() %%%

	iter2 = 0;

	while s_gap > ep / 3
		%%% Shift() %%%
		z = z + s_gap / log(m);
		% z = z - s_gap / (10 * sqrt(M));
			%%% NewtonStep() %%%
			s = c - [A' * y(1 : n) + y(n + 1 : end); y(n + 1 : end); -y(n + 1 : end); y(1 : n); -y(1 : n)];
				%%% Solve(A*S_inv2*A', -As_inv) %%%
				s_inv = s.^(-1);
				r = [-A * s_inv(1 : m) - s_inv(3 * m + 1 : 3 * m + n) + s_inv(3 * m + n + 1 : end);
					 -s_inv(1 : m) - s_inv(m + 1 : 2 * m) + s_inv(2 * m + 1 : 3 * m)];
				r = r + (M / (b' * y - z)) * b;

				s_inv2 = s.^(-2);
				D1 = spdiags(s_inv2(1 : m), 0, m, m);
				D2 = spdiags(s_inv2(3 * m + 1 : 3 * m + n) + s_inv2(3 * m + n + 1 : M), 0, n, n);
				d13 = s_inv2(1 : m) + s_inv2(m + 1 : 2 * m) + s_inv2(2 * m + 1 : 3 * m);
				d13 = spdiags(d13.^(-1), 0, m, m);

				Ms = A * D1 * A' + D2 - A * D1 * d13 * D1 * A';
				Q = A * D1 * d13;
				v = sqrt(M) / (b' * y - z) * b;

				% x1 = cmgSolver(Ms, r(1 : n) - Q * r(n + 1 : end));
				[x1, garbage] = pcg(Ms, r(1 : n) - Q * r(n + 1 : end), ep1);
				x1 = [x1; d13 * (r(n + 1 : end) - D1 * A' * x1)];

				% x2 = cmgSolver(Ms, v(1 : n) - Q * v(n + 1 : end));
				[x2, garbage] = pcg(Ms, v(1 : n) - Q * v(n + 1 : end), ep1);
				x2 = [x2; d13 * (v(n + 1 : end) - D1 * A' * x2)];

				x = x1 - (x2 * (x2' * r)) / (1 + v' * x2);
				%%% End Solve() %%%
			y = y + (1 - 1 / (20 * (sqrt(M) + 1))) * x;
			%%% End NewtonStep() %%%
		s_gap = b' * y - z;
		%%% End Shift() %%%
		iter2 = iter2 + 1;
	end
	iterations_of_shift = iter2

	s = c - [A' * y(1 : n) + y(n + 1 : end); y(n + 1 : end); -y(n + 1 : end); y(1 : n); -y(1 : n)];
	%%% Solve(A*S_inv2*A', As_inv) %%%
	s_inv = s.^(-1);
	r = [A * s_inv(1 : m) + s_inv(3 * m + 1 : 3 * m + n) - s_inv(3 * m + n + 1 : end);
		 s_inv(1 : m) + s_inv(m + 1 : 2 * m) - s_inv(2 * m + 1 : 3 * m)];
	r = r - (M / s_gap) * b;

	s_inv2 = s.^(-2);
	D1 = spdiags(s_inv2(1 : m), 0, m, m);
	D2 = spdiags(s_inv2(3 * m + 1 : 3 * m + n) + s_inv2(3 * m + n + 1 : M), 0, n, n);
	d13 = s_inv2(1 : m) + s_inv2(m + 1 : 2 * m) + s_inv2(2 * m + 1 : 3 * m);
	d13 = spdiags(d13.^(-1), 0, m, m);

	Ms = A * D1 * A' + D2 - A * D1 * d13 * D1 * A';
	Q = A * D1 * d13;
	v = sqrt(M) / s_gap * b;

	% x1 = cmgSolver(Ms, r(1 : n) - Q * r(n + 1 : end));
	[x1, garbage] = pcg(Ms, r(1 : n) - Q * r(n + 1 : end), ep1);
	x1 = [x1; d13 * (r(n + 1 : end) - D1 * A' * x1)];

	% x2 = cmgSolver(Ms, v(1 : n) - Q * v(n + 1 : end));
	[x2, garbage] = pcg(Ms, v(1 : n) - Q * v(n + 1 : end), ep1);
	x2 = [x2; d13 * (v(n + 1 : end) - D1 * A' * x2)];

	v = x1 - (x2 * (x2' * r)) / (1 + v' * x2);

	s1 = spdiags(s_inv2(1 : m), 0, m, m);
	v1 = v(1 : n);
	v2 = v(n + 1 : end);
	%%% end Solve() %%%

	x = [s1 * A' * v1 + s1 * v2;
		 spdiags(s_inv2(m + 1 : 2 * m), 0, m, m) * v2;
		 -spdiags(s_inv2(2 * m + 1 : 3 * m), 0, m, m) * v2;
		 spdiags(s_inv2(3 * m + 1 : 3 * m + n), 0, n, n) * v1;
		 -spdiags(s_inv2(3 * m + n + 1 : end), 0, n, n) * v1];
	x = s_inv - x;
	x = x / M;
	x = x / ((1 / s_gap) + (b' * v / (s_gap^(2))));

	x = round(x(m, 1));
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
