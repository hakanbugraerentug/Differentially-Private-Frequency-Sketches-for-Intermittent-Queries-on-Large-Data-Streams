function f_est = DP_count_query(x_q, hash_params, eps_DP)

global C U;

% sketch type:
% 1: count sketch
% 2: winterfell (take the minimum used only)
% 3: winterfallen (exclude the maximum used only, unless every cell is used same many times)
% 4: catsandogs: uses all the cells
% 5: 

% answer type:
% 1: median
% 2: mean


d = size(hash_params{1}, 1);
w = hash_params{1}(1, end);

% initialise the set

n_q = length(x_q);

f_est = zeros(2, n_q);

% first step: find the cells that are to be added noise
E = zeros(d, n_q);
E_cell = cell(1, d);
for i = 1:d
    for k = 1:n_q
        x = x_q(k);
        E(i, k) = my_hash(x, hash_params{1}(i, :));
    end
    E_cell{i} = unique(E(i, :));
end

% second step: add noise to the cells
for i = 1:d
    for j = E_cell{i}
        noise = laprnd(d/eps_DP);
        C(i, j) = C(i, j) + noise;
        U(i, j) = U(i, j) + 1;
    end
end

% third step: add noise to the cells and calculate median
for k = 1:n_q
    x = x_q(k);
    w = zeros(1, d);
    C_for_median = zeros(1, d);
    for i = 1:d
        j = my_hash(x, hash_params{1}(i, :));
        s = (-1)^my_hash(x, hash_params{2}(i, :));
        C_for_median(i) = s*C(i, j);
        w(i) = U(i, j);
    end
    
    % first estimate: simply the median
    f_est(1, k) = median(C_for_median);
    
    % second estimate: weight w.r.t. noise
    v = 1./w;
    sum_v = sum(v);
    [C_sorted, b] = sort(C_for_median);
    
    c = 0; temp_sum = 0; j = 1;
    while c == 0
        temp_sum = temp_sum + v(b(j));
        if temp_sum > sum_v/2
            c = 1;
        else
            j = j + 1;
        end
    end
    f_est(2, k) = C_sorted(j);
end


