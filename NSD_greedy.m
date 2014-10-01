function [aM, aG, sim_sparsity, dt] = NSD_greedy(A, B, preiters, iters, alpha)
% alignment matching aM
% alignment graph aG
% struct with timings dt

dt = [];
dt.components = 0.0;
% prepare uniform vectors as initial conditions
% and compute for preiters steps (alpha = 1.0)
t0 = clock;
n = size(A, 1);
m = size(B, 1);
vecA = ones(n, 1) ./ n;
vecB = ones(m, 1) ./ m;

vecsA = inpowers(A, preiters, vecA);
vecsB = inpowers(B, preiters, vecB);

ivecA = vecsA{preiters + 1};
ivecB = vecsB{preiters + 1};

ivecsA = {};
ivecsB = {};

ivecsA{1} = ivecA;
ivecsB{1} = ivecB;
dt_preiters = etime(clock, t0);  
dt.preiters = dt_preiters;

% compute similarity matrix 
[X, dt_iters] = NSD_rank_compute(A, B, alpha, iters, ivecsA, ivecsB);
dt.iters = dt_iters;

% compute matching
[M, dt_match] = greedy_match(X);
dt.match = dt_match;

% compute alignment graph (and the corresponding matching)
[aM, aG, dt_align] = align(A, B, M);
dt.align = dt_align;

sim_sparsity = 1.0;
end

