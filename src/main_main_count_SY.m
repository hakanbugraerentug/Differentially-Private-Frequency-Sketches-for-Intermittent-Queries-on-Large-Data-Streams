M = 1000000; % universe size
T = 1000000; % stream size
L0 = 100000; % query size
% number of repetitions
Rep = 3;
eps_DP = 1; % epsilon for DP
eps_DP0_vec = [inf];
L_eps = length(eps_DP0_vec);

alpha_vec = [1.6];
L_a = length(alpha_vec);

w_vec = [2000];
d_vec = [5];
w_by_wu = 10;
R = 10;
L_w = length(w_vec);
L_d = length(d_vec);
alg_set = [5];

query_type_vec = [1 2];
L_q = length(query_type_vec);

for i1 = 1:L_eps
    eps_DP0 = eps_DP0_vec(i1);
    for i2 = 1:L_a
        alpha = alpha_vec(i2);
        for i3 = 1:L_w
            w = w_vec(i3);
            w_u = floor(w/w_by_wu);
            for i4 = 1:L_d
                d = d_vec(i4);
                for i5 = 1:L_q
                    query_type = query_type_vec(i5);
                    main_count_as_function(M, T, L0, Rep, query_type, ...
                        alpha, eps_DP, eps_DP0, d, w, w_u, R, alg_set);
                end
            end
        end
    end
end
