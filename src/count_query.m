function f_est = count_query(x, hash_params, sketch_type, eps_DP)

global C U V;

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

% initialise the set
C_for_median = zeros(1, d);
weight = ones(1, d);

switch sketch_type
    
    case 1 % count sketch
        for i = 1:d
            j = my_hash(x, hash_params{1}(i, :));
            s = (-1)^my_hash(x, hash_params{2}(i, :));
            C_for_median(i) = s*C(i, j);
        end
    case 2 % winterfell
        % determine the minimum:
        u_min = inf;
        index_vec = zeros(1, d);
        u_vec = zeros(1, d);
        for i = 1:d
            j = my_hash(x, hash_params{1}(i, :));
            index_vec(i) = j;
            u_vec(i) = U(i, j);
            u_min = min(u_min, U(i, j));
        end
        
        % update and use the cells that are used for the minimum number of
        % times
        k = 0;
        for i = 1:d
            if u_vec(i) == u_min
                j = index_vec(i);
                s = (-1)^my_hash(x, hash_params{2}(i, :));
                k = k + 1;
                noise = laprnd(d/eps_DP);
                C_for_median(k) = s*C(i, j) + noise;
                
                % update the tables
                U(i, j) = U(i, j) + 1;
                C(i, j) = C(i, j) + noise;
                weight(k) = U(i, j);
            end
        end
        C_for_median = C_for_median(1:k);
        weight = weight(1:k);
        
    case 3 % winterfallen
        % determine the minimum and the maximum:
        u_min = inf;
        u_max = 0;
        index_vec = zeros(1, d);
        u_vec = zeros(1, d);
        for i = 1:d
            j = my_hash(x, hash_params{1}(i, :));
            index_vec(i) = j;
            u_vec(i) = U(i, j);
            u_min = min(u_min, U(i, j));
            u_max = max(u_max, U(i, j));
        end
        
        k = 0;
        for i = 1:d
            if u_vec(i) < u_max || u_max == u_min
                j = index_vec(i);
                s = (-1)^my_hash(x, hash_params{2}(i, :));
                k = k + 1;
                noise = laprnd(d/eps_DP);
                C_for_median(k) = s*C(i, j) + noise;
                
                % update the tables
                U(i, j) = U(i, j) + 1;
                C(i, j) = C(i, j) + noise;
                weight(k) = U(i, j);
            end
        end
        
        C_for_median = C_for_median(1:k);
        weight = weight(1:k);
        
    case 4 % catsanddogs
        for i = 1:d
            j = my_hash(x, hash_params{1}(i, :));
            s = (-1)^my_hash(x, hash_params{2}(i, :));
            noise = laprnd(d/eps_DP);
            C_for_median(i) = s*C(i, j) + noise;
            
            % update the tables
            C(i, j) = C(i, j) + noise;
            U(i, j) = U(i, j) + 1;
            weight(i) = U(i, j);
        end
        
    case 5 % avoidcollusions
        collusion = zeros(1, d);
        i_vec_consider = zeros(1, d);
        for i = 1:d
            j = my_hash(x, hash_params{1}(i, :));
            s = (-1)^my_hash(x, hash_params{2}(i, :));
            noise = laprnd(d/eps_DP);
            C_for_median(i) = s*C(i, j) + noise;
            C(i, j) = C(i, j) + noise;
            
            j_h = my_hash(j, hash_params{3}(i, :));
            collusion(i) = U(i, j_h);
            
            if V(i, j_h) == j
                i_vec_consider(i) = 1;
            end
        end
        if sum(i_vec_consider) > 0
            max_collusion = max(collusion(i_vec_consider == 1));
            min_collusion = min(collusion(i_vec_consider == 1));
            diff_collusion = max_collusion - min_collusion;
            selection_vec = (max_collusion - collusion) >= 0.01*diff_collusion...
                | (max_collusion*ones(1, d) <= 1000)...
                | (i_vec_consider == 0);
        else
            selection_vec = 1:d;
        end
        C_for_median = C_for_median(selection_vec);
        weight = weight(selection_vec);
        
    case 6 % winterwillfall
        
        % determine the minimum and the maximum:
        u_max = 0;
        index_vec = zeros(1, d);
        u_vec = zeros(1, d);
        for i = 1:d
            j = my_hash(x, hash_params{1}(i, :));
            index_vec(i) = j;
            u_vec(i) = U(i, j);
            u_max = max(u_max, U(i, j));
        end
        
        k = 0;
        max_discarded = 0;
        for i = 1:d            
            if u_vec(i) < u_max || max_discarded == 1
                j = index_vec(i);
                s = (-1)^my_hash(x, hash_params{2}(i, :));
                k = k + 1;
                noise = laprnd(d/eps_DP);
                C_for_median(k) = s*C(i, j) + noise;
                
                % update the tables
                U(i, j) = U(i, j) + 1;
                C(i, j) = C(i, j) + noise;
                weight(k) = U(i, j);
            else
                max_discarded = 1;
            end
        end
        
        C_for_median = C_for_median(1:k);
        weight = weight(1:k);
        
end

% return the answer (either by taking the median or mean)
f_est(1) = median(C_for_median);
f_est(2) = mean(C_for_median);


% w_eight = (1./(1+weight))/sum(1./(1+weight));

v = sum(1./weight);
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
f_est(3) = C_sorted(j);

% f_est(3) = sum(C_for_median.*w_eight);

if sketch_type == 5
    for i = 1:d
        j = my_hash(x, hash_params{1}(i, :));
        j_h =  my_hash(j, hash_params{3}(i, :));
        [U(i, j_h), b_temp] = max([U(i, j_h), f_est(1)]);
        if b_temp == 2
            V(i, j_h) = j;
        end
    end
end

