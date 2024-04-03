function [soln_numeric, soln_analytic] = winner_combin_rows(u0, u_vec, epsilon, M)

% This function numerically compares all subsets of rows given the number 
% of times a laplace(d/epsilon) is added to the value of the cell.
d = length(u_vec);
[u_sorted, ~] = sort(u_vec);
X = zeros(d, M);
X_median = zeros(1, d);

for j = 1:d
    X(j, :) = sum((-1).^(rand(u_sorted(j), M) < 0.5).*exprnd(d/epsilon, u_sorted(j), M), 1);
end
for j = 1:d
    X_median(j) = mean(abs(mean(X(1:j, :), 1)));
    % X_median(j) = std(mean(X(1:j, :), 1));
end
[~, temp] = min(X_median);

soln_numeric = u_sorted(1:temp);

% analytical solution
temp_vec = cumsum((d/epsilon)*u_sorted)./((1:d).^2);

[~, temp] = min(temp_vec);
soln_analytic = u_sorted(1:temp);
