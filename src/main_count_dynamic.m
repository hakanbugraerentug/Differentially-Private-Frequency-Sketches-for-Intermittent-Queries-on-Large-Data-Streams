M = 1000000; % universe size
T = 1000000; % stream size
L = 1000000; % query size
% number of repetitions
eps_DP = 1; % epsilon for DP

alpha_vec = [1.6];
L_a = length(alpha_vec);

w_vec = [2000];
d_vec = [5];
R = 1;
L_w = length(w_vec);
L_d = length(d_vec);

P_vec = [10000];
L_P = length(P_vec);

query_type_vec = [2];
L_q = length(query_type_vec);

for i1 = 1:L_a
    alpha = alpha_vec(i1);
    for i2 = 1:L_w
        w = w_vec(i2);
        for i3 = 1:L_d
            d = d_vec(i3);
            for i4 = 1:L_q
                for i5 = 1:L_P
                    P = P_vec(i5);
                    query_type = query_type_vec(i4);
                    [mae{i1, i2, i3, i4, i5}, ...
                        mpe1{i1, i2, i3, i4, i5}, ...
                        mpe2{i1, i2, i3, i4, i5}]...
                        = main_count_dynamic_as_function(M, T, L, alpha, query_type, eps_DP, d, w, P, R);
                end
            end
        end
    end
end

