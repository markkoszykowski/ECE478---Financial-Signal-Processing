function [S0, V0] = wrapper(V, S0, u, d, r, K, n, m, tolerance, p)
% Computes value of underlying stock and asset at an initial provided time
% through a Monte Carlo simulation
    % if probability of H coin flip (up stock movement) provided generate 
    % random walks with that, otherwise use risk neutral probabilities
    if ~exist("p", "var")
        [Sn, ~] = generatemrandompaths(S0, u, d, ...
            ((1 + r) - d) / (u - d), (u - (1 + r)) / (u - d), n, m, tolerance);
    else
        [Sn, ~] = generatemrandompaths(S0, u, d, p, 1 - p, n, m, tolerance);
    end
    Vn = V(Sn, K);

    % 'bring back' values of underlying stock and asset to initial time
    S0 = mean(Sn(:, end)) / (1 + r) ^ n;
    V0 = mean(Vn(:, end)) / (1 + r) ^ n;
end