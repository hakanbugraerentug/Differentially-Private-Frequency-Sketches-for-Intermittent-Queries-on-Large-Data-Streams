% M = 1000;
% N = 100000;
% k = 100;
% y = zeros(1, M);
% for m = 1:M
%     v = rand(1, N);
%     v_sort = sort(v, 'ascend');
%     y(m) = (k-1)/v_sort(k);
% end
% 
% z = zeros(1, M);
% 
% for m = 1:M
%     z_temp = 0;
%     
%     min_vec = ones(1, k);
%     max_min_vec = max(min_vec);
%     for i = 1:N
%         u = rand;
%         if u < max_min_vec
%             z_temp = z_temp + 1/max_min_vec;
%             min_vec_temp = sort(min_vec, 'ascend');
%             min_vec = [min_vec_temp(1:k-1) u];
%             max_min_vec = max(min_vec);
%         end
%     end
%     z(m) = z_temp;
% end
% MAE_1 = mean(abs(y - N));
% MAE_2 = mean(abs(z - N));
% 
% MSE_1 = mean((y - N).^2/N^2);
% MSE_2 = mean((z - N).^2/N^2);
% 
% %%
% figure;
% subplot(1, 2, 1);
% histogram(y, 50);
% title(sprintf('k-min value, MSE = %.3f', MAE_1));
% subplot(1, 2, 2);
% histogram(z, 50);
% title(sprintf('k-min with geometric increments, MSE = %.3f', MAE_2));
% sgtitle('Histograms of 1000 estimates of the true cardinality 10^5');
% 
% 
% %%
% T = 1000; k = 100; t_vec = (k+1):T;
% V1 = [zeros(1,k) k^2*((t_vec-k+1).*t_vec)/((k-1)^2*(k-2))];
% V2 = zeros(1, T);
% for t = k+1:T
%     tau_vec = k+1:t;
%     V2(t) = sum(tau_vec/(k-1) - 1);
% end
% figure;
% subplot(1, 2, 1);
% plot(1:T, V1, 1:T, V2);
% legend('standard', 'alternative');
% 
% subplot(1, 2, 2);
% plot(V1./V2);
% title('V1/V2');
% hold on;
% plot(exp(1)*ones(1, T), 'r');
% hold off;
% 
%% hyperloglog?
B = 64;
A_exp = 10;
A = 2^A_exp;
alpha = 0.697;

M = 1000;
N = 10000;

MemCost1 = A*ceil(log2(B-A_exp));

y = zeros(1, M);

for m = 1:M
    if mod(m, 100) == 0
        disp(m);
    end
    C = zeros(1, A);
    for i = 1:N
        c = randi([1, A]);
        u = randi([0 1], 1, B-A_exp);
        k = find(u, 1);
        if isempty(k)
            k = B-A_exp+1;
        end
        if k > C(c)
            C(c) = k;
        end 
    end
    
    y(m) = alpha*A*1/mean(2.^(-C));
end

Cost = A*log(B - A_exp);

%
A_exp = 10;
A2 = 2^A_exp;

z = zeros(1, M);
z2 = zeros(1, M);

MemCost2 = 32 + A2*ceil(log2(B-A_exp));

for m = 1:M
    if mod(m, 100) == 0
        disp(m);
    end
    
    C = zeros(1, A2);
    prob_change = 1;
    z_temp = 0;    
    for i = 1:N
        c = randi([1, A2]);
        u = randi([0 1], 1, B-A_exp);
        k = find(u, 1);
        if isempty(k)
            k = B-A_exp+1;
        end
        
        if k > C(c)
            z_temp = z_temp + 1/prob_change;
            prob_change = prob_change - 2.^(-C(c))/A2;
            C(c) = k;
            prob_change = prob_change + 2.^(-C(c))/A2;
        end 
    end  
    z(m) = z_temp;
    z2(m) = alpha*A2*1/mean(2.^(-C));
end

%%
histogram(y, 50, 'Normalization', 'probability', 'FaceColor', 'g');
[H, Bins] = hist(y, 50);

hold on;
histogram(z, Bins, 'Normalization', 'probability', 'FaceColor','k');
hold off;
legend('HLL', 'new');

MSE_1 = mean((y - N).^2/N^2);
MSE_2 = mean((z - N).^2/N^2);
MSE_3 = mean((z2 - N).^2/N^2);
MSE_4 = mean(((z2+2*z)/3- N).^2/N^2);

MAE_1 = mean(abs(y - N)/N);
MAE_2 = mean(abs(z - N)/N);
MAE_3 = mean(abs(z2 - N)/N);
MAE_4 = mean((abs(z2+2*z)/3 - N)/N);