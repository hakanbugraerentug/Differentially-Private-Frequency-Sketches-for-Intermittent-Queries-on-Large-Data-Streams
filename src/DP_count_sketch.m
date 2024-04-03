function [outputs] = DP_count_sketch(X, X_q, d, w, hash_params, sketch_type, eps_DP, eps_DP0)

global C;
global U;
global V;

C = zeros(d, w);
U = zeros(d, w);
if sketch_type == 5
    w_u = hash_params{3}(1, 4);
    U = zeros(d, w_u);
    V = zeros(d, w_u);
end

% stream size
T = length(X);

% Query size
L = length(X_q);

% construct the count table
for t = 1:T
    x = X(t);
    update_C(x, hash_params); 
end

% add noise for DP algorithms
if sketch_type > 1
    C = C + laprnd(d/eps_DP0, d, w);
end

% answer queries
f_est = zeros(L, 3);
for i = 1:L
    x = X_q(i);
    f_est(i, :) = count_query(x, hash_params, sketch_type, eps_DP);
end

outputs.f_est = f_est;
outputs.U = U;
outputs.C = C;

