I = 3;
K = 10000;
y = zeros(I, K);

w = 1*(1:I).^2;
v = 1./w;
sum_v = sum(v);
for i = 1:I
    y(i, :) = laprnd(w(i), 1, K);
end


x_est_0 = mean(y);
x_est_1 = median(y);

x_est_2 = zeros(1, K);
for k = 1:K
    [y_sorted, b] = sort(y(:, k));
    c = 0; temp_sum = 0; j = 1;
    while c == 0
        temp_sum = temp_sum + v(b(j));
        if temp_sum > sum_v/2
            c = 1;
        else
            j = j + 1;
        end
    end
    x_est_2(k) = y_sorted(j);
end

plot(x_est_1);
hold on;
plot(x_est_2, 'r');
hold off;

E0 = mean(abs(x_est_0));
E1 = mean(abs(x_est_1));
E2 = mean(abs(x_est_2));

disp([E0, E1, E2]);