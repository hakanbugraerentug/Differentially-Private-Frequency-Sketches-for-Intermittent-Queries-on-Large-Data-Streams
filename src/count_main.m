clear all; clc; close all;

n = 100; M = 10000; epsilon = 0.1;
X_mean = zeros(1, n);
X_median = zeros(1, n);
for i = 1:n
    X = (-1).^(rand(i, M) < 0.5).*exprnd(i/epsilon, i, M);
    % X = i*randn(i, M);
    X_mean(i) = mean(abs(mean(X)));
    X_median(i) = mean(abs(median(X)));
end
plot(X_mean)
hold on;
plot(X_median);
hold off;

%%

clear all; clc; close all;
d = 5;
u_max = 6;
winner = zeros(u_max, u_max, u_max, u_max, u_max);
X_median = zeros(1, d);
M = 100000; epsilon = 0.1;
winner_rows = cell(1, u_max^d);
winner_vec = zeros(1, u_max^d);
u_comb = zeros(u_max^d, d);
i = 0;
for i1 = 1:u_max
    for i2 = i1:u_max
        for i3 = i2:u_max
            for i4 = i3:u_max
                for i5 = i4:u_max
                    i = i + 1;
                    i_vec = [i1, i2, i3, i4, i5];
                    u_comb(i, :) = i_vec;
                    temp_rows = winner_combin_rows(i_vec, 1, M);
                    winner_rows{i} = temp_rows;
                    winner_vec(i) = length(temp_rows);
                end
            end
        end
    end
end

                    