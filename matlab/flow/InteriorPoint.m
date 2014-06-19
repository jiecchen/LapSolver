function x = InteriorPoint(A, b, c, lmin, T, y, ep)
	n = length(b);
	m = length(c);

	tic;
	[y, s, s_gap, z] = FindCentralPath(A, b, c, lmin, T, y);
	toc;
	fprintf('FindCentralPath\n');

	tic;
	while s_gap > (ep / 3)
		[y, s, s_gap, z] = Shift(A, b, c, y, s_gap, z);
	end
	toc;
	fprintf('Shift\n');

	A_new = [A, -b * ones(1, m)];
	s_new = [s; s_gap * ones(m, 1)];
	s_min = min(s_new);

	max_A = max(max(abs(A)));
	max_b = max(abs(b));
	max_c = max(abs(c));
	U = max([max_A, max_b, max_c]);

	%%%
	%% ep4 = min(1, (s_min * sqrt(m)) / (T * U * n));
	%%%

	S_new_inv = sparse(diag(s_new.^(-1)));
	S_new_inv2 = sparse(diag(s_new.^(-2)));

	%%%
	%% v = Solve (A_new * S_new_inv2 * A_new') X = A_new * S_new_inv * ones(2*m, 1)
	%%%
	tic;
	v = amdSolver(A_new * S_new_inv2 * A_new', A_new * S_new_inv * ones(2*m, 1));
	toc;
	fprintf('amdSolver\n');


	S_inv = sparse(diag(s.^(-1)));
	S_inv2 = sparse(diag(s.^(-2)));

	x = (S_inv * ones(m, 1) - S_inv2 * A' * v) / ((m / s_gap) + (m / s_gap^2) * b' * v);
end

function [y, s, s_gap, z] = Shift(A, b, c, y, s_gap, z)
	m = length(c);

	z = z + (s_gap / (10 * sqrt(m)));
	A_new = [A, -b * ones(1, m)];
	c_new = [c; -z * ones(m, 1)];

	y = NewtonStep(A_new, c_new, y);
	s = c - A' * y;
	s_gap = b' * y - z;
end

function y = NewtonStep(A, c, y)
	m = length(c);

	s = c - A' * y;
	S_inv = sparse(diag(s.^(-1)));
	S_inv2 = sparse(diag(s.^(-2)));

	%%%
	%% dy = Solve (A * S_inv2 * A') X = -A * S_inv * ones(m, 1)
	%%%
	dy = amdSolver(A * S_inv2 * A', -A * S_inv * ones(m, 1));

	y = y + dy * (1 - 1 / (20 * (sqrt(m) + 1)));
end

function [y, s, s_gap, z] = FindCentralPath(A, b, c, lmin, T, y)
	m = length(c);

	s = c - A' * y;
	S_inv = sparse(diag(s.^(-1)));
	s_gap = m;

	b_new = A * S_inv * ones(m, 1);
	z = b_new' * y - s_gap;

	%% what should the bound be??? 
	bound =  40 * T * m * norm(b_new) / sqrt(lmin);

	while s_gap < bound
		[y, s, s_gap, z] = Unshift(A, b_new, c, y, s_gap, z);
	end

	z = b' * y - bound * norm(b) / norm(b_new);
	s = c - A' * y;
	s_gap = b' * y - z;
end

function [y, s, s_gap, z] = Unshift(A, b_new, c, y, s_gap, z)
	m = length(c);

	z = z - s_gap / (10 * sqrt(m)); % changed from sqrt(m) to log(m)
	A_new = [A, -b_new * ones(1, m)];
	c_new = [c; -z * ones(m, 1)];

	y = NewtonStep(A_new, c_new, y);

	s = c - A' * y;
	s_gap = b_new' * y - z;
end





















