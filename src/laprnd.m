function x = laprnd(b, m, n)

if nargin == 1
    m = 1;
    n = 1;
elseif nargin == 2
    m = 1;
    n = m;
end

x = exprnd(b, m, n).*(-1).^(rand(m, n) < 0.5);