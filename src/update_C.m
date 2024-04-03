function update_C(x, hash_params)

global C;

d = size(hash_params{1}, 1);
for i = 1:d
    j = my_hash(x, hash_params{1}(i, :));
    s = (-1)^my_hash(x, hash_params{2}(i, :));
    C(i, j) = C(i, j) + s;
end