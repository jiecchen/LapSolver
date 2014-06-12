function x = InteriorPoint(A, b, c, lmin, T, y, ep)
% function x = InteriorPoint(A, b, c, lmin, T, y, ep)
% A is an nxm matrix with linearly independent rows
% b is a length n vector
% c is a length m vector
% lmin > 0 is lower bound on the eigenvalues of AA'
% T > 0 is an upper bound on the absolute values of coordinates of the feasible dual polytope
% y is a length n vector satisfying A'y < c
% 0 < ep < 1

	n = length(b);
	m = length(c);

	[y, s, sg, z] = FindCentralPath(A, b, c, lmin, T, y);
	s = c - A' * y;
	sg = b' * y - z;

	while sg > (ep / 3)
		[y, s, sg, z] = Shift(A, b, c, y, sg, z);
	end

	An = [A, -b * ones(1, m)];
	sn = [s; sg * ones(m, 1)];
	smin = min(sn);
	% T = max(max(y), max(s));
	U = max([max(max(abs(A))), max(abs(b)), max(abs(c))]);
	ep4 = min(1, (smin / (T * U)) * (sqrt(m) / n));

	S = diag(s);
	Sn = diag(sn);

	% v = pcg(An * inv(Sn^2) * An', An * inv(Sn) * ones(2 * m, 1), ep4);
	v = inv(An * inv(Sn^2) * An')* An * inv(Sn) * ones(2 * m, 1);

	x = (inv(S) * ones(m, 1) - inv(S^2) * A' * v) / ((m / sg) + (m / (sg^2)) * b' * v);
end

function [y, s, sg, z] = Shift(A, b, c, y, sg, z)

	m = length(c);

	z = z + sg / (10 * sqrt(m));
	An = [A, -b * ones(1, m)];
	cn = [c; -z * ones(m, 1)];

	y = NewtonStep(An, cn, y);
	s = c - A' * y;
	sg = b' * y - z;
end

function y = NewtonStep(A, c, y) 

	m = length(c);

	s = c - A' * y;
	S = diag(s);
	dy = inv(A * inv(S^2) * A') * -A * inv(S) * ones(m, 1);
	y = y + dy * (1 - 1 / (20 * (sqrt(m) + 1)));
end

function [y, s, sg, z] = FindCentralPath(A, b, c, lmin, T, y)

	m = length(c);

	s = c - A' * y;
	S = diag(s);
	sg = m;
	bn = A * inv(S) * ones(m, 1);
	z = bn' * y - sg;

	% T = max(max(y), max(s));
	bound = (40 / sqrt(lmin)) * T * m * norm(bn);

	while sg < bound
		[y, s, sg, z] = Unshift(A, b, c, y, sg, z);
	end

	z = b' * y - (40 / sqrt(lmin)) * T * m * norm(b);
	sg = b' * y - z;
end

function [y, s, sg, z] = Unshift(A, b, c, y, sg, z)

	m = length(c);

	z = z - sg / (10 * sqrt(m));
	An = [A, -b * ones(1, m)];
	cn = [c; -z * ones(m, 1)];

	y = NewtonStep(An, cn, y);

	s = c - A'* y;
	sg = b' * y - z;
end













