function emission_map = generate_pseudo_emissions(size_x,size_y,alpha)
% Generate emission map
gaus = @(x,mu,sig)exp(-(((x-mu).^2)/(2*sig.^2)));
gaus2D = @(x, y, mu,sig)kron(gaus(x, mu, sig), gaus(y, mu, sig));

x_vec = linspace(-2, 2, size_x);
y_vec = linspace(-2, 2, size_y);
dist = alpha * gaus2D(x_vec, y_vec, 0, 1);

emission_map = gamrnd(dist, 1);
emission_map = bsxfun(@rdivide, emission_map, sum(emission_map,2));
emission_map = 100 * 1000 * emission_map; % let the total emissions be 100
emission_map = reshape(emission_map, size_x, size_y);
end

