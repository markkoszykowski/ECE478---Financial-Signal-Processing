function [dW, W] = wiener(N, delta, M, includewiener)
% Generates 'M' Wiener discretized processes of length 'N' with time
% interval delta
    dW = normrnd(0, sqrt(delta), [M, N]);
    if exist("includewiener", "var")
        W = zeros(M, N+1);
        for n = 1:N
            W(:, n+1) = sum(dW(:, 1:n), 2);
        end
    end
end