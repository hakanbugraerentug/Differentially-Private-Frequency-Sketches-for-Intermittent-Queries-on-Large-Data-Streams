clear all; clc; close all; fc = 0;


U = 10^6;
N = 10^6;

epsilon = 1;

n_vec = 1000;
L_n = length(n_vec);

Num_of_trials = 10000;

k_vec = [10 20 50 100];
L_k = length(k_vec);

n_est =cell(L_n, L_k);
e_r = cell(L_n, L_k);

for i1 = 1:L_n
    n = n_vec(i1);
    for i2 = 1:L_k
        k = k_vec(i2);
        disp([n k]);
        % determine m
        m_min = ceil(k/epsilon - 1);
        m = m_min;
        c = k*log(m+1) - epsilon*(m+1-k) < -log(U);
        while c == 0
            m = m+1;
            c = k*log(m+1) - epsilon*(m+1-k) < -log(U);
        end
        
        % initialize the hash table
        S = N/n;
        e_r_temp = zeros(Num_of_trials, S);
        for trial = 1:Num_of_trials
            y_ext = sort([ones(1, k) rand(1, m)], 'ascend');
            y = y_ext(1:k);
            for s = 1:S
                y_ext = sort([y rand(1, n)], 'ascend');
                y = y_ext(1:k);
                n_est_temp = (k-1)/y(end) - s*m;
                e_r_temp(trial, s) = abs(n_est_temp - n*s)/(n*s);
                
                % add artificial hashes
                y_ext = sort([y rand(1, m)]);
                y = y_ext(1:k);
            end
        end
        e_r{i1, i2} = mean(e_r_temp);
        
    end
end

%%
fc = fc + 1; figure(fc);
legend_cell = cell(1, L_k);
for i2 = 1:L_k
    legend_cell{i2} = sprintf('k = %d', k_vec(i2));
end

for i1 = 1:L_n
    n = n_vec(i1);
    
    subplot(1, L_n, i1);
    S = N/n_vec(i1);
    for i2 = 1:L_k-1
        plot(1:S, e_r{i1, i2});
        hold on;
    end
    hold off;
    legend(legend_cell(1:3));
    xlabel('query number');
    title(sprintf('query at every %d unique members', n));
    ylabel('relative error');
end
