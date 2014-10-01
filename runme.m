
data = load('data/dmela.mat'); A = data.A;
data = load('data/scere.mat'); B = data.A;
preiters = 2; iters = 10; alpha = 0.8;
[M, G, ~, ~] = NSD_greedy(A, B, preiters, iters, alpha);
