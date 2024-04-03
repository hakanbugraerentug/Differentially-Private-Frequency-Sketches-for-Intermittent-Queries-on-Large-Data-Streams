M = 1000;
N = 100000;
k = 100;
y = zeros(1, M);
for m = 1:M
    v = rand(1, N);
    v_sort = sort(v, 'ascend');
    y(m) = (k-1)/v_sort(k);
end

z = zeros(1, M);
for m = 1:M
    z_temp = 0;
    min_vec = ones(1, k);
    max_min_vec = max(min_vec);
    for i = 1:N
        u = rand;
        if u < max_min_vec
            min_vec_temp = sort(min_vec, 'ascend');
            min_vec = [min_vec_temp(1:k-1) u];
            max_min_vec = max(min_vec);
            z_temp = z_temp + 1/max_min_vec;
        end
    end
    z(m) = z_temp;
end

MAE_1 = mean(abs(y - N));
MAE_2 = mean(abs(z - N));

figure;
subplot(1, 2, 1);
histogram(y, 50);
title(sprintf('k-min value, MAE = %.3f', MAE_1));
subplot(1, 2, 2);
histogram(z, 50);
title(sprintf('k-min with geometric increments, MAE = %.3f', MAE_2));
sgtitle('Histograms of 1000 estimates of the true cardinality 10^5');
