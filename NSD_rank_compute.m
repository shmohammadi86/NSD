function [X, dt] = NSD_rank_compute(A, B, alpha, iters, vecsA, vecsB, sigmas)
% NSD_rank_compute Computes the similarity of two undirected graphs using the
% decomposition approach.
% 
% Input arguments:
% - A, B: the adjacency matrices of the two graphs.
% - alpha: alpha parameter of the IsoRank algorithm.
% - iters: the number of iterations.
% - vecsA, vecsB: cell arrays containing the vectors which will participate 
%     in composing the preferences matrix H. 
%     Summing outer products vecsB{i} * vecsA{i}' is meant to be 
%     a decomposition of the preferences matrix H, when no sigmas vector is given
%     (typically from NMF or some specific user constructed decomposition).
%     Summing sigmas(i) * vecsB{i} * vecsA{i}'  is meant to be 
%     a decomposition of H when sigmas are given, typically in the SVD case.
%  - sigmas: a vector of sigma values (typically a number of the leading
%     singular values coming from SVD). Optional parameter.
%
% Output arguments:
% - X: the matrix with the similarity scores (similarity matrix).
%     Note that element X(i,j) is the similarity score of node i in B 
%     and node j in A; if B has m nodes and A has n nodes then X is an 
%     m x n matrix.
% - dt: the time in seconds for the operation.


% Giorgos Kollias and Shahin Mohammadi
% Department of Computer Science, Purdue University

A = max(A, A');
B = max(B, B');


s = size(vecsA, 2);

ivecsA = {};
ivecsB = {};

m = size(B, 1);
n = size(A, 1);
X = zeros(m, n);

t0 = clock;    
if  ~exist('sigmas', 'var') || isempty(sigmas)
    for ss=1:s
       ivecsA{ss} = vecsA{ss};
       ivecsB{ss} = vecsB{ss};
    end
else
    for ss=1:s
       ivecsA{ss} = sigmas(ss) * vecsA{ss};
       ivecsB{ss} = vecsB{ss};
    end
end  

k = iters + 1;


for ss=1:s
    vecsA = inpowers(A, iters, ivecsA{:, ss});
    vecsB = inpowers(B, iters, ivecsB{:, ss});
    factor = 1.0 - alpha;
    for i = 1:k-1
        term = factor * vecsB{i} * vecsA{i}';
        factor = factor * alpha;
        X = X + term;
    end
    factor = alpha ^ iters;
    term = factor * vecsB{k} * vecsA{k}';
    X = X + term;

end

dt = etime(clock, t0);


end

