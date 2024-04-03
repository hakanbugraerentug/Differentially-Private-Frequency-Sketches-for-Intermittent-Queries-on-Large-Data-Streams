function [mae, mpe1, mpe2] = main_count_dynamic_as_function(M, T, L, alpha, query_type, eps_DP, d, w, P, R)

close all;

% generate data and query set
X = zipf_rand(M, alpha, T);
if query_type == 1
    X_q = zipf_rand(M, 1, L);
elseif query_type == 2
    X_q = randperm(L);
elseif query_type == 3
    X_q = 1:L;
end

% select a prime number between M and 2M
p_vec = primes(2*M);
p_vec = p_vec(p_vec > M);
p = p_vec(ceil(length(p_vec)*rand));

outputs = cell(5, R);
mae = cell(5, 3);
mpe1 = cell(5, 3);
mpe2 = cell(5, 3);

%% 
for r = 1:R
    disp(r);
    % randomly generate the hash functions
    for i = 1:d
        a = ceil(rand*p);
        b = ceil(rand*p);
        hash_params{1}(i, :) = [a, b, p, w];
        a = ceil(rand*p);
        b = ceil(rand*p);
        hash_params{2}(i, :) = [a, b, p, 2];
    end
    % Run the algorithms
    outputs{r} = DP_count_dynamic(X, X_q, M, P, d, w, hash_params, eps_DP);
    F_true = outputs{r}.F_true;
    F_est = outputs{r}.F_est;
    
    mae{1}(r, :) = abs(F_est(1, :) - F_true);
    mpe1{1}(r, :) = abs(F_est(1, :) - F_true)./max(1, F_true);
    mpe2{1}(r, :) = abs(F_est(1, :) - F_true)./max(1, max(F_est(1, :), F_true));
        
    mae{2}(r, :) = abs(F_est(2, :) - F_true);
    mpe1{2}(r, :) = abs(F_est(2, :) - F_true)./max(1, F_true);
    mpe2{2}(r, :) = abs(F_est(2, :) - F_true)./max(1, max(F_est(2, :), F_true));
end
