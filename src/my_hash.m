function y = my_hash(x, hash_param)

a = hash_param(1);
b = hash_param(2);
p = hash_param(3);
w = hash_param(4);

y = mod(mod(a*x + b, p), w) + 1;