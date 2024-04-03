mu0 = 0.5:0.01:1;
k = 1;
mu = mu0*k;
delta = 1 - k./(2*mu);

X = (exp(-delta)./((1-delta).^(1-delta))).^mu;
Y = exp(-(delta.^2).*mu/2);

plot(mu0, X); hold on; plot(mu0, Y, 'r'); hold off;