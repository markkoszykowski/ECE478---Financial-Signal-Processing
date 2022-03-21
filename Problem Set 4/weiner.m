function [W] = weiner(N, delta, M)
    W = zeros(M, N+1);
    for n = 1:N
        W(:, n+1) = W(:, n) + normrnd(0, sqrt(delta), [M, 1]);
    end
end