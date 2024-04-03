function main_count_as_function(M, T, L0, Rep, query_type, alpha, eps_DP, eps_DP0, d, w, w_u, R, alg_set)

close all;

L = Rep*L0;
sort_opt = 1;

% generate data and query set
X = zipf_rand(M, alpha, T);
if query_type == 1
    X_q = zipf_rand(M, 1, L0);
elseif query_type == 2
    X_q = randperm(L0);
elseif query_type == 3
    X_q = 1:L0;
end

[~, ind_sort] = sort(X_q);

% true frequencies
f_true = histc(X, 1:M)';
f_true_q = f_true(X_q);

X_q = repmat(X_q, 1, Rep);

% select a prime number between M and 2M
p_vec = primes(2*M);
p_vec = p_vec(p_vec > M);
p = p_vec(ceil(length(p_vec)*rand));

outputs = cell(5, R);
cum_mae = cell(5, 3);
cum_mpe1 = cell(5, 3);
cum_mpe2 = cell(5, 3);
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
        hash_params{3}(i, :) = [1, -1, w_u, w_u];
    end
    % Run the algorithms
    for i = alg_set
        hash_params_input = hash_params;
        w_input = w;
        if i == 1 || i == 4 % for count sketch or cats and dogs
            hash_params_input{1}(:, 4) = w + w_u;
            w_input = w + w_u;
        end
        [outputs{i, r}] = DP_count_sketch(X, X_q, d, w_input, ...
            hash_params_input, i, eps_DP, eps_DP0);
    end
end

%%
for r = 1:R
    for i = alg_set
        for j = 1:3
            f_est_temp = max(0, outputs{i, r}.f_est(:, j));
            if sort_opt == 1
                ind_vec_rep = reshape(L0*repmat(0:(Rep-1), L0, 1), L, 1)' ...
                    + repmat(ind_sort, 1, Rep);
                f_est_temp = f_est_temp(ind_vec_rep);
                f_true_q_rep = repmat(f_true_q(ind_sort), Rep, 1);
            else
                f_true_q_rep = repmat(f_true_q, Rep, 1);
            end
            
            mae{i, j}(r, :) = abs(f_est_temp - f_true_q_rep);
            mpe1{i, j}(r, :) = abs(f_est_temp - f_true_q_rep)./max(1, f_true_q_rep);
            mpe2{i, j}(r, :) = abs(f_est_temp - f_true_q_rep)./max(1, max(f_est_temp, f_true_q_rep));

            cum_mae{i, j}(r, :) = cumsum(mae{i, j}(r, :));
            cum_mpe1{i, j}(r, :) = cumsum(mpe1{i, j}(r, :));
            cum_mpe2{i, j}(r, :) = cumsum(mpe2{i, j}(r, :));
        end
    end
end

%
legends1 = {'CS', 'WF', 'WFn', 'CD', 'AC'};
legends2 = {'first pass', 'second pass', 'third pass', 'fourth pass'};

figure(1);
subplot(2, 2, 1);
hold on;
for i = alg_set
    plot(mean(cum_mae{i, 1}, 1));
end
hold off;
xlabel('query no'); title('cumulative MAE');
grid on;

x = mean(cum_mae{1, 1}, 1) - mean(cum_mae{5, 1}, 1);
subplot(2, 2, 3);
plot(x(1:L0));
hold on;
for i = 2:Rep
    plot(x((i-1)*L0+1:i*L0) - x((i-1)*L0));
end
% plot(x(2*L0+1:3*L0) - x(2*L0), 'k');
hold off;
xlabel('query no'); title('difference in cumul. MAE');
grid on;

subplot(2, 2, 2);
hold on;
for i = alg_set
    plot(mean(cum_mpe1{i, 1}, 1));
end
hold off;
xlabel('query no'); title('cumulative MPE');
grid on;
legend(legends1(alg_set), 'Location','southeast');

x = mean(cum_mpe1{1, 1}, 1) - mean(cum_mpe1{5, 1}, 1);
subplot(2, 2, 4);
plot(x(1:L0));
hold on;
for i = 2:Rep
    plot(x((i-1)*L0+1:i*L0) - x((i-1)*L0));
end
hold off;
xlabel('query no'); title('difference in cumul. MAE');
grid on;
legend(legends2(1:Rep), 'Location','southeast');
t = sgtitle(sprintf('M=%d, T=%d, w=%d, wu = %d, d=%d, alpha=%02d, epsDP0=%03d, query type = %d', M, T, w, w_u, d, 10*alpha, 100*eps_DP0, query_type));
t.FontSize = 8;
filename_to_save = sprintf('fig_M_%d_T_%d_query_type_%d_alpha_%02d_eps_DP_%03d_eps_DP0_%03d_d_%d_w_%d_wu_%d_R_%d',...
    M, T, query_type, 10*alpha, 100*eps_DP, 100*eps_DP0, d, w, w_u, R);
fname = strcat(filename_to_save, '.eps');
print(fname, '-depsc');
% save(filename_to_save);
close all;
