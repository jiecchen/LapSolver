% function [t1, t2, siz] = test()
function test(n)

	A = grid2(n, n);
	x = maxflow(A, 1, n*n, ones(size(a2u(A), 2), 1))
 
end

% 	n = 4;
% 	iters = 40;

% 	t1 = zeros(iters, 1);
% 	t2 = zeros(iters, 1);
% 	siz = zeros(iters, 1);

% 	for i = 1:iters

% 		A = delGraph(n);
% 		% c = randi([1, 10], size(a2u(A), 2), 1);
% 		c = ones(size(a2u(A), 2), 1);
% 		tic;
% 		x = maxflow1(A, 1, n, c);
% 		t1(i, 1) = toc;
% 		tic;
% 		y = graphmaxflow(triu(A), 1, n, 'capacity', c); 
% 		t2(i, 1) = toc;

% 		siz(i, 1) = length(A);

% 		n = n + 100

% 	end
% end

% function [t3, siz2] = test()

% 	n = 4004
% 	iters = 200
% 	t3 = zeros(iters, 1);
% 	siz2 = zeros(iters, 1);

% 	for i = 1:iters

% 		A = delGraph(n);
% 		c = ones(size(a2u(A), 2), 1);
% 		tic;
% 		y = graphmaxflow(triu(A), 1, n, 'capacity', c); 
% 		t3(i, 1) = toc;

% 		siz2(i, 1) = length(A);

% 		n = n + 1000
% 	end


% 	% if x == y
% 	% 	same = true
% 	% else
% 	% 	same = false
% 	% end
% end