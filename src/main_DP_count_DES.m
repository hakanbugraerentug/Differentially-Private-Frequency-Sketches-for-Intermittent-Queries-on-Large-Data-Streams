% discrete event simulation for the use and keep algorithm
clear all; clc; close all; fc = 0;
% main code for comparing the count sketches
M = 1000000; % universe size
T = 1000000; % stream size
N = 100000; % query range

alpha = 1.6;

% Parameters of the experiments
eps_DP_vec = [1]; L_e = length(eps_DP_vec); % epsilon for DP
d_vec = [3 5 9]; L_d = length(d_vec); % number of rows
n_p_vec = [100 100000]; L_n = length(n_p_vec); % periods for query
prop_const = 100;
n_Q_vec = n_p_vec/prop_const; % query size vec;
T_min = 100000; % minimum stream size after which the queries arrive
NQ = (T-T_min)/prop_const; % total number of queries
u_vec = [1 100]; L_u = length(u_vec); % user type vector, this is the percentage of the top N entries that are queried

w = 2000; % number of columns
R = 10; % number of experiments per each

%% For hashing: select a prime number between M and 2M
p_vec = primes(2*M);
p_vec = p_vec(p_vec > M);
p = p_vec(ceil(length(p_vec)*rand));

%% generate data and query set
X = zipf_rand(M, alpha, T);

%%

Error = cell(L_e, L_d, L_n, L_u, 2);

for i1 = 1:L_e    
    eps_DP = eps_DP_vec(i1); % get the epsilon parameter
    for i2 = 1:L_d
        d = d_vec(i2); % get d
        for i3 = 1:L_n
            n_p = n_p_vec(i3); % get the period
            n_Q = n_Q_vec(i3); % get the query size
            
            for i4 = 1:L_u
                u = u_vec(i4);

                % initialize the error vectors
                error_1 = zeros(R, NQ);
                error_2 = zeros(R, NQ);
                error_3 = zeros(R, NQ);

                for r = 1:R
                    disp([i1, i2, i3, i4, r]);
                    % initialise the true frequency vector
                    f_true = zeros(1, T);
                    % randomly generate the hash functions
                    for i5 = 1:d
                        % for the location in the hash table
                        a = ceil(rand*p);
                        b = ceil(rand*p);
                        hash_params{1}(i5, :) = [a, b, p, w];
                        % for the sign of the hashed value
                        a = ceil(rand*p);
                        b = ceil(rand*p);
                        hash_params{2}(i5, :) = [a, b, p, 2];
                    end
                    
                    % initialise the hash table
                    C = zeros(d, w);
                    q = 0; % initialize the query number
                    for t = 1:T
                        % get the next element and update the true freq
                        x = X(t);
                        f_true(x) = f_true(x) + 1;
                        % update the hash table
                        for i = 1:d
                            j = my_hash(x, hash_params{1}(i, :));
                            s = (-1)^my_hash(x, hash_params{2}(i, :));
                            C(i, j) = C(i, j) + s;        
                        end
                        
                        if mod(t, n_p) == 0 && t > T_min
                            % generate the query set
                            Q = randsample(N*u/100, n_Q);
                            
                            % answer the queries
                            E = repmat({zeros(1, n_Q)}, 1, d);
                            c = zeros(1, d);
                            
                            for i6 = 1:n_Q
                                q = q + 1; % increment the query
                                x_q = Q(i6); % get the element to be queried
                                med_vec = zeros(1, d);
                                for i7 = 1:d
                                    j = my_hash(x_q, hash_params{1}(i7, :));
                                    s = (-1)^my_hash(x_q, hash_params{2}(i7, :));
                                    if sum(E{i7} == j) == 0 % if this is a new cell
                                        c(i7) = c(i7) + 1; 
                                        E{i7}(c(i7)) = j; % add the column index
                                        C(i7, j) = C(i7, j) + laprnd(d/eps_DP); % add noise
                                    end
                                    % add the value in the vector for
                                    % median calculation
                                    med_vec(i7) = s*C(i7, j);
                                end
                                % median estimate
                                f_x_est = max(0, median(med_vec));
                                % calculate the error
                                f_x_true = f_true(x_q);
                                error_1(r, q) = abs(f_x_est-f_x_true);
                                error_2(r, q) = abs(f_x_est-f_x_true)/max(f_x_true, 1);
                            end
                        end % if query time
                    end % t 
                end % r
                Error{i1, i2, i3, i4, 1} = mean(error_1);
                Error{i1, i2, i3, i4, 2} = mean(error_2);
            end % i4
        end % i3
    end % i2
end % i1


%%
legend_cell = cell(1, L_d);
for i = 1:L_d
    legend_cell{i} = sprintf('d = %d', d_vec(i));
end
ylabel_string = {'cumul. mean absolute error', 'cumul. mean relative error'};

for error_type = 1:2
    for i1 = 1:L_e
        eps_DP = eps_DP_vec(i1);
        fc = fc + 1; figure(fc);
        for i3 = 1:L_n
            n_p = n_p_vec(i3);
            for i4  = 1:L_u
                u = u_vec(i4);
                subplot(L_n, L_u, (i3-1)*L_u +i4);
                for i2 = 1:L_d
                    d = d_vec(i2);
                    plot(cumsum(Error{i1, i2, i3, i4, error_type}));
                    hold on;
                end
                hold off;
                if i3 == L_n && i4 == L_u 
                    legend(legend_cell);
                end
                if i3 == L_n
                    xlabel('number of queries');
                end
                if i4 == 1
                    ylabel(ylabel_string(error_type));
                end
                title(sprintf('n_p = %d, u = %.2f', n_p, u*N/T));
            end
        end
    end
end