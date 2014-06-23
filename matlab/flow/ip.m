function x = ip(A, b, c, y, T, lmin, ep)
	% function x = ip(A, b, c, y, T, lmin, ep)

	m = length(c);

	%%%%%%%%%%%%%% FIND CENTRAL PATH
	s = c - A' * y;
	s_gap = m;
	S_inv = sparse(diag(s.^(-1)));
	b_ = A * S_inv * ones(m, 1);
	z = b_' * y - s_gap;

	bound = 40 * T * m * norm(b_) / sqrt(lmin);

	iter1 = 0;

	while s_gap < bound
		%z = z - s_gap / (10 * sqrt(m)); % begin unshift
		z = z - s_gap / log(m);
		tmpA = [A, -b_ * ones(1, m)];
		tmpc = [c; -z * ones(m, 1)];

		s = tmpc - tmpA' * y;
		tmp1 = tmpA * sparse(diag(s.^(-2))) * tmpA';
		tmp2 = -tmpA * sparse(diag(s.^(-1))) * ones(2 * m, 1);
		%%%%%%%%%%%%%% SOLVE(tmp1, tmp2)
		dy = tmp1 \ tmp2;
		%%%%%%%%%%%%%%
		y = y + dy * (1 - 1 / (20 * (sqrt(m) + 1)));

		s = c - A' * y;
		s_gap = b_' * y - z; % end unshift


		iter1 = iter1 + 1;


	end
	z = b' * y - 40 * T * m * norm(b) / sqrt(lmin);
	s_gap = b' * y - z;
	%%%%%%%%%%%%%%

	iter2 = 0;
	while s_gap > ep / 3
		%z = z + s_gap / (10 * sqrt(m)); % begin shift
		z = z + s_gap / log(m);
		tmpA = [A, -b * ones(1, m)];
		tmpc = [c; -z * ones(m, 1)];

		s = tmpc - tmpA' * y;
		tmp1 = tmpA * sparse(diag(s.^(-2))) * tmpA';
		tmp2 = -tmpA * sparse(diag(s.^(-1))) * ones(2 * m, 1);
		%%%%%%%%%%%%%% SOLVE(tmp1, tmp2)
		dy = tmp1 \ tmp2;
		%%%%%%%%%%%%%%
		y = y + dy * (1 - 1 / (20 * (sqrt(m) + 1)));

		s = c - A' * y;
		s_gap = b' * y - z; % end shift


		iter2 = iter2 + 1;

	end
	a_new = [A, -b * ones(1, m)];
	s_new = [s; s_gap * ones(m, 1)];
	tmp1 = a_new * sparse(diag(s_new.^(-2))) * a_new';
	tmp2 = a_new * sparse(diag(s_new.^(-1))) * ones(2 * m, 1);

	%%%%%%%%%%%%%% SOLVE(tmp1, tmp2)
	v = tmp1 \ tmp2;
	%%%%%%%%%%%%%%

	tmp1 = sparse(diag(s.^(-1))) * ones(m, 1) - sparse(diag(s.^(-2))) * A' * v;
	tmp2 = (1 / s_gap) + (b' * v / (s_gap^(2)));
	x = tmp1 / (m * tmp2);

	iterations_of_unshift = iter1
	iterations_of_shift = iter2

end


% function y = NewtonStep(A, c, y, m)
% 	s = c - A' * y;
% 	tmp1 = A * sparse(diag(s.^(-2))) * A';
% 	tmp2 = -A * sparse(diag(s.^(-1))) * ones(2 * m, 1);

% 	%%%%%%%%%%%%%% SOLVE(tmp1, tmp2)
% 	dy = tmp1 \ tmp2;
% 	%%%%%%%%%%%%%%

% 	y = y + dy * (1 - 1 / (20 * (sqrt(m) + 1)));
% end













