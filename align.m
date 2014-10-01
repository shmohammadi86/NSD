function [aM, aG, dt] = align(A, B, M)

t0 = clock;

% ensure adjacency, symmetric
A(A > 0) = 1.0;
A = max(A, A');

B(B > 0) = 1.0;
B = max(B, B');

% we assume m <=n 
n = size(A, 1);
m = size(B, 1);
[mb, ma] = find(M);

[s_mb, ind] = sort(mb); % sort matching nodes of smaller graph (B) 
p_ma = ma(ind); % permute matching nodes of larger graph (A) accordingly
C = A(p_ma, p_ma);
D = B + C;
p = spdiags(D, 0);
D = D - spdiags(p, 0, m, m);
[ii, jj] = find(D == 2); % B and the permuted A have elements at these coords

aM = sparse(s_mb, p_ma, 1.0, m, n); % matching with the new numbering
aG = sparse(ii, jj, 1.0, m, m); % alignment graph

dt = etime(clock, t0);

end

