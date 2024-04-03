clear all; clc; close all; fc = 0;

% main code for comparing the count sketches
M = 1000000; % universe size
T = 1000000; % stream size
L0 = 100000; % query size
% number of repetations
Rep = 1;
L = Rep*L0;
eps_DP = 1; % epsilon for DP
eps_DP0 = inf; % epsilon for DP
sort_opt = 1;
alpha = 1.6; % parameter of the zipf distribution

% generate data and query set
X = zipf_rand(M, alpha, T);

X_q = zipf_rand(M, 1, L0);
% X_q = sort(X_q);
% X_q = [1:L];
% X_q = randperm(L0);

[X_sort, ind_sort] = sort(X_q);

% true frequencies
f_true = histc(X, 1:M)';
f_true_q = f_true(X_q);

X_q = repmat(X_q, 1, Rep);

d = 3;
w = 1000;

% select a prime number between M and 2M
p_vec = primes(2*M);
p_vec = p_vec(p_vec > M);
p = p_vec(ceil(length(p_vec)*rand));

w_u = floor(w/inf);

% Number of MC runs
R = 10;

alg_set = [1 4 6];

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

%%
legends = {'CS', 'WF', 'WFn', 'CD', 'AC'};
fc = fc + 1; figure(fc);
subplot(2, 2, 1);
hold on;
for i = alg_set
    plot(mean(cum_mae{i, 1}, 1));
end
hold off;
xlabel('query no'); title('cumulative MAE');
grid on;
legend('CS', 'CD', 'WW');
%
x = mean(cum_mae{1, 1}, 1) - mean(cum_mae{6, 1}, 1);
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

x = mean(cum_mpe1{1, 1}, 1) - mean(cum_mpe1{6, 1}, 1);
subplot(2, 2, 4);
plot(x(1:L0));
hold on;
for i = 2:Rep
    plot(x((i-1)*L0+1:i*L0) - x((i-1)*L0));
end
hold off;
xlabel('query no'); title('difference in cumul. MAE');
grid on;
legend('first pass', 'second pass', 'third pass');

% subplot(2, 3, 3);
% hold on;
% for i = alg_set
%     plot(mean(cum_mpe2{i, 1}));
% end
% hold off;
% xlabel('query no'); title('cumulative MPE2');
% grid on;
% legend(legends{alg_set});
% 
% x = mean(cum_mpe2{1, 1}) - mean(cum_mpe2{5, 1});
% subplot(2, 3, 6);
% plot(x(1:L0));
% hold on;
% for i = 2:Rep
%     plot(x((i-1)*L0+1:i*L0) - x((i-1)L0));
% end
% hold off;
% xlabel('query no'); title('difference in cumul. MAE');
% grid on;


%%

% fc = fc + 1; figure(fc);
% subplot(3, 3, 1);
% hold on;
% for i = alg_set
%     plot(mean(cum_mae{i, 1}));
% end
% hold off;
% xlabel('query no'); title('cumulative MAE');
% subplot(3, 3, 2);
% hold on;
% for i = alg_set
%     plot(mean(cum_mpe1{i, 1}));
% end
% hold off;
% xlabel('query no'); title('cumulative MPE');
% subplot(3, 3, 3);
% hold on;
% for i = 1:4
%     plot(mean(cum_mpe2{i, 1}));
% end
% hold off;
% xlabel('query no'); title('cumulative MPE2');
% 
% subplot(3, 3, 4);
% hold on;
% for i = alg_set
%     plot(mean(cum_mae{i, 2}));
% end
% ylabel('mean of the values');
% hold off;
% xlabel('query no'); title('cumulative MAE');
% subplot(3, 3, 5);
% hold on;
% for i = alg_set
%     plot(mean(cum_mpe1{i, 2}));
% end
% hold off;
% xlabel('query no'); title('cumulative MPE');
% subplot(3, 3, 6);
% hold on;
% for i = alg_set
%     plot(mean(cum_mpe2{i, 2}));
% end
% hold off;
% xlabel('query no'); title('cumulative MPE2');
% 
% subplot(3, 3, 7);
% hold on;
% for i = alg_set
%     plot(mean(cum_mae{i, 3}));
% end
% ylabel('mean of the values');
% hold off;
% xlabel('query no'); title('cumulative MAE');
% subplot(3, 3, 8);
% hold on;
% for i = alg_set
%     plot(mean(cum_mpe1{i, 3}));
% end
% hold off;
% xlabel('query no'); title('cumulative MPE');
% subplot(3, 3, 9);
% hold on;
% for i = alg_set
%     plot(mean(cum_mpe2{i, 3}));
% end
% hold off;
% xlabel('query no'); title('cumulative MPE2');
% legend('CS', 'WF', 'WFn', 'CD');