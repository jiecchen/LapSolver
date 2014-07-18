function test(n)

	A = grid2(n, n);
	tic;
	x = maxflow(A, 1, n*n, ones(size(a2u(A), 2), 1))
	toc;
	tic;
	x1 = maxflow1(A, 1, n*n, ones(size(a2u(A), 2), 1))
	toc;
 
end
