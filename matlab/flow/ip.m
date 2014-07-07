function x = ip(A, b, c, y, T, lmin, ep)
	% function x = ip(A, b, c, y, T, lmin, ep)

	% n = size(A, 1);
	m = size(A, 2);

	s = c - A' * y;
	s_gap = m;
	s_inv = s.^(-1);
	b_ = A * s_inv;
	z = b_' * y - s_gap;

	bound = 40 * T * m * norm(b_) / sqrt(lmin);

	iter1 = 0;

	% tmpA = [A, -b_ * sparse(ones(1, m))];
	bbt = b_ * b_';
	while s_gap < bound
		z = z - s_gap / log(m);
		%z = z - s_gap / (10 * sqrt(m));
		% tmpc = [c; -z * sparse(ones(m, 1))];
		% s = tmpc - tmpA' * y;
		s1 = c - A' * y;
		s2 = -z + b_' * y;
		s2 = s2 * ones(m, 1);
		s = [s1; s2];

		s_inv = s.^(-1);
		s_inv2 = s_inv.^(2);

		thinga = spdiags(s_inv2(1:m), 0, m, m);
		thingb = sum(s_inv2(m + 1:end));
		tmp1 = A * thinga * A' + thingb * bbt;

		thinga = s_inv(1:m);
		thingb = sum(s_inv(m + 1:end));
		tmp2 = -A * thinga + thingb * b_;

		[dy, garbage] = pcg(tmp1, tmp2, 1 / (20 * (sqrt(m) + 1)));

		dy = sparse(dy);
		y = y + dy * (1 - 1 / (20 * (sqrt(m) + 1)));

		s = c - A' * y;
		s_gap = b_' * y - z;

		iter1 = iter1 + 1;



	end
	z = b' * y - 40 * T * m * norm(b) / sqrt(lmin);
	s_gap = b' * y - z;

	iter2 = 0;

	% tmpA = [A, -b * sparse(ones(1, m))];
	bbt = b * b';

	while s_gap > ep / 3
		z = z + s_gap / log(m);
		%z = z + s_gap / (10 * sqrt(m));

		% tmpc = [c; -z * sparse(ones(m, 1))];

		% s = tmpc - tmpA' * y;
		s1 = c - A' * y;
		s2 = -z + b' * y;
		s2 = s2 * ones(m, 1);
		s = [s1; s2];

		s_inv = s.^(-1);
		s_inv2 = s_inv.^(2);

		thinga = spdiags(s_inv2(1:m), 0, m, m);
		thingb = sum(s_inv2(m + 1:end));
		tmp1 = A * thinga * A' + thingb * bbt;

		thinga = s_inv(1:m);
		thingb = sum(s_inv(m + 1:end));
		tmp2 = -A * thinga + thingb * b;

		[dy, garbage] = pcg(tmp1, tmp2, 1 / (20 * (sqrt(m) + 1)));
		dy = sparse(dy);
		y = y + dy * (1 - 1 / (20 * (sqrt(m) + 1)));

		s = c - A' * y;
		s_gap = b' * y - z;

		iter2 = iter2 + 1;

	end

	% a_new = [A, -b * sparse(ones(1, m))];
	s_new = [s; s_gap * sparse(ones(m, 1))];

	s_inv = s_new.^(-1);
	s_inv2 = s_inv.^(2);

	thinga = spdiags(s_inv2(1:m), 0, m, m);
	thingb = sum(s_inv2(m + 1:end));
	tmp1 = A * thinga * A' + thingb * bbt;

	thinga = s_inv(1:m);
	thingb = sum(s_inv(m + 1:end));
	tmp2 = -A * thinga + thingb * b;

	[v, garbage] = pcg(tmp1, tmp2, 0.2);
	v = sparse(v);

	s_inv = s.^(-1);
	s_inv2 = s_inv.^(2);

	tmp1 = s.^(-1) - spdiags(s_inv2(:), 0, m, m) * A' * v;
	tmp2 = (1 / s_gap) + (b' * v / (s_gap^(2)));

	x = tmp1 / (m * tmp2);

	iterations_of_unshift = iter1
	iterations_of_shift = iter2
end