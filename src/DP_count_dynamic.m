function [outputs] = DP_count_dynamic(X, X_q, M, P, d, w, hash_params, eps_DP)

global C;
global U;

% construct the count table
C = zeros(d, w);
U = zeros(d, w);

% stream size
T = length(X);
NQ = length(X_q);

% calculate the query set size per query
nQ = NQ/P;

% calculate the query period
pQ = T/P;

F_est = zeros(2, nQ);

% period counter
p = 0;

% true counts
f_true = zeros(1, M);
F_true = zeros(1, NQ);

for t = 1:T
    if mod(t, 100000) == 0
        disp(t);
    end
    x = X(t);
    % update the true array
    f_true(x) = f_true(x) + 1;
    
    update_C(x, hash_params);
    % check for queries
    if mod(t, pQ) == 0
        p = p + 1;
        % get the query set
        q_ind = (p-1)*nQ+1: p*nQ;
        x_q = X_q(q_ind);
        % F_est is a matrix whose rows contain responses using different
        % methods
        F_est(:, q_ind) = DP_count_query(x_q, hash_params, eps_DP);
        F_true(q_ind) = f_true(x_q);
    end
end

outputs.F_est = F_est;
outputs.F_true = F_true;
outputs.U = U;
outputs.C = C;