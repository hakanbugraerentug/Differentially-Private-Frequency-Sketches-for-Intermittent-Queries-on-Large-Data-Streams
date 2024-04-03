% clc; clear all; close all;
% D = 20;
% d_vec = 2*(1:D) - 1;
% 
% x = zeros(1, D);
% y = zeros(1, D);
% 
% u_vec = [0 100 1000];
% L_u = length(u_vec);
% M = 10000;
% for j = 1:L_u
%     u = u_vec(j);
%     for i = 1:D
%         d = d_vec(i);
%         B = [zeros(floor(d/2+1), 1); u*ones(d-floor(d/2+1), 1)];
%         x(i) = mean(abs(median(repmat(B, 1, M) + laprnd(d, d, M), 1) ));
%         y(i) = mean(abs(laprnd(d, 1, M)));
%     end
%     plot(d_vec, x);
%     hold on;
% end
% plot(d_vec, d_vec, '-.');
% hold off;
% grid on;
% legend('u = 0', 'u = 1', 'u = 10', 'mae of noisy median');
% ylabel('mae of median');
% xlabel('d');
% title('mae of median of noisy added vector values');
% 
% %%
% clc; clear all; % close all;
% d_vec = 2*(1:20) - 1;
% d_vec = 1:50;
% D = length(d_vec);
% epsilon = 1;
% w = 1000;
% f_2 = 1000;
% lambda_target = 1/3;
% kappa_vec = zeros(1, D);
% 
% for i = 1:D
%     d = d_vec(i);
%     u_vec = ones(1, d);
%     kappa_min_sq = 2*max(u_vec.*(d^2/epsilon^2) + f_2^2/w);
%     kappa_sq = kappa_min_sq;
%     if -d/2 + log(2)*(d/2) > log(lambda_target)
%         % if -d/8 > log(lambda_target)
%         kappa_vec(i) = inf;
%     else
%         
%         c = 0;
%         while c == 0
%             lambda = (f_2^2/(w) + d^2/(epsilon^2)*mean(u_vec))/kappa_sq;
%             lambda_med = -d/2*((1 - 2*lambda) - log(2*(1 - lambda)));
%             % lambda_med = -d/8*((1 - 2*lambda)^2/(1 - lambda));
%             
%             if lambda_med < log(lambda_target)
%                 c = 1;
%             else
%                 kappa_sq = kappa_sq + 1;
%             end
%             
%         end
%         kappa_vec(i) = sqrt(kappa_sq);
%     end
% end
% hold on;
% plot(d_vec, kappa_vec);
% grid on;
% 

%
clear all; clc; close all;
nQ_vec = 1:2000;
w_vec = [1000 5000 10000];
L_n = length(nQ_vec);
L_w  = length(w_vec);

EmQ1 = zeros(L_w, L_n);
EmQ2 = zeros(L_w, L_n);
for i = 1:L_w
    w = w_vec(i);
    for j = 1:L_n
        nQ = nQ_vec(j);
        x = 1:nQ;
        EmQ1(i, j) = sum(min(1, w*betainc(1/w,x,nQ-x+1)));
        temp_alpha = w*betainc(1/w,x-1,nQ-x+2);
        EmQ2(i, j) = sum(min(1, max(0, temp_alpha - temp_alpha.^2*(w-1)/(2*w))));
    end
    
end

colorlist = {'k', 'r', 'b'};
for i = 1:L_w
    % plot(nQ_vec, EmQ1(i, :), nQ_vec, EmQ2(i, :), colorlist{i});
    plot(nQ_vec, EmQ1(i, :), colorlist{i});
    hold on;
end
hold off;
legend('w = 1000', 'w = 10000');
title('upper bound for expected value of $\hat{m}_Q/d$', 'Interpreter', 'Latex');
xlabel('$n_{Q}$', 'Interpreter', 'Latex'); ylabel('$E(\hat{m}_{Q}/d)$', 'Interpreter', 'Latex');
grid on;

%% exact probabilities

% This section calculates the exact probabilities 
% P_b(n, m): The probability of having cell size at most b out of n
% multinomial samples into m buckets with equal probability.
%
% source for the recursion can be found at
% http://www.nebrija.es/~jmaestro/esa/papers/JDA2011.pdf

% clear all;
% clc; % close all;

N = 2000;
W = 10000;
temp_prev_vec = zeros(N, N);

f_vec = [0 cumsum(log(1:N))];

m_vec = [1000 5000 10000];
L_m = length(m_vec);
EmQ = zeros(L_m, N);

m_target = m_vec(1);
c = 1;
for m = 1:W
    disp(m);
    prob_mtx = ones(N, N);
    if m == 1
        for b = 1:N
            for n = 1:N
                prob_mtx(n, b) = (b>=n);
            end
        end
    else
        for b = 1:N
            
            % temp_vec = ones(N, 1);
            % temp_vec = ones(N, 1);
            
 %           tic;
            n_vec = (b+1):N;
            temp_factor_vec = [1 prev_prob_mtx(n_vec(2:end)-1-b, b)'];
            
            sub_vec = f_vec(n_vec) - f_vec(b+1) - f_vec(n_vec-b)...
                        + log(temp_factor_vec) + (n_vec-1-b)*log(m-1) - (n_vec-1)*log(m);
            
            max_sub_vec =  max(sub_vec);
            cum_sub_vec = log(cumsum(exp(sub_vec - max_sub_vec))) + max_sub_vec;
            
            prob_mtx(n_vec, b) = max(0, 1 - exp(cum_sub_vec));
            
%                toc;     
%             tic;
%             for n = (b+1):min(N, m*b)
%                 % if b < n
%                     if n > b+1
%                         temp_factor = prev_prob_mtx(n-1-b, b);
%                     else
%                         temp_factor = 1;
%                     end
%                     log_minus_term = f_vec(n-1+1) - f_vec(b+1) - f_vec(n-1-b+1)...
%                         + log(temp_factor) + (n-1-b)*log(m-1) - (n-1)*log(m);
%                     
%                     [temp_vec(n), max_ind] = max([0, temp_vec(n-1) - exp(log_minus_term)]);
%                     
%                 % else
%                 %     temp_vec(n) = 1; 
%                 % end  
%                 
%                
%                 if max_ind == 1
%                     disp([m, b, n]);
%                     break;
%                 end
%             end
%            toc;
 %           disp(sum(abs(temp_vec(n_vec) - temp_vec2(n_vec))));
        end
    end
    % disp(m);
    prev_prob_mtx = prob_mtx;
    
    if m == m_target
        EmQ(c, :) = 1+sum(1 - prob_mtx, 2); 
        c = c + 1;
        if c <= L_m
            m_target = m_vec(c);
        end
    end
end

%%
for i = 1:L_w
    plot(1:N, EmQ(i, :), colorlist{i});
    hold on;
end
hold off;
legend('w = 1000', 'w = 5000', 'w = 10000');
title('expected value of $\hat{m}_Q/d$', 'Interpreter', 'Latex');
xlabel('$n_{Q}$', 'Interpreter', 'Latex'); ylabel('$E(\hat{m}_{Q}/d)$', 'Interpreter', 'Latex');
grid on;
