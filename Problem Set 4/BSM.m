function [c] = BSM(tau, x, K, r, sigma)
    d_plus = (1 / (sigma * sqrt(tau))) * (log(x / K) + (r + (sigma ^ 2 / 2)) * tau);
    d_minus = (1 / (sigma * sqrt(tau))) * (log(x / K) + (r - (sigma ^ 2 / 2)) * tau);
    c = x .* normcdf(d_plus) - K * exp(-r * tau) * normcdf(d_minus);
end

