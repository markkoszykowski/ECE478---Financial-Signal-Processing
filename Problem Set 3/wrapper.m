function [S0, V0] = wrapper(V, S0, u, d, r, K, n, m, tolerance, p)
    if ~exist("p", "var")
        [Sn, ~] = generatemrandompaths(S0, u, d, ...
            ((1 + r) - d) / (u - d), (u - (1 + r)) / (u - d), n, m, tolerance);
    else
        [Sn, ~] = generatemrandompaths(S0, u, d, p, 1 - p, n, m, tolerance);
    end
    Vn = V(Sn, K);

    S0 = mean(Sn(:, end)) / (1 + r) ^ n;
    V0 = mean(Vn(:, end)) / (1 + r) ^ n;
end