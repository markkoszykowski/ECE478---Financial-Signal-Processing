function [X, dW, paths] = dXanddiscard(X, n, beta, gamma, dt, dW, paths, verbose)
% Determines next value of artbitrary asset X using a discretized
% stochastic differential equation with time coefficient beta and Wiener
% process coefficient gamma
    dX = beta .* dt + gamma .* dW(:, n-1);
    Xn_1 = dX + X(:, n-1);

    % ensure that value of asset is never less than or equal to zero
    X(:, n) = Xn_1 .* (Xn_1 > 0);

    % discard impractical paths
    X = X((Xn_1 > 0), :);
    dW = dW((Xn_1 > 0), :);
    % print an exception for discarded paths
    if verbose
        if any(Xn_1 <= 0, "all")
            for m = find(Xn_1 <= 0).'
                disp("Path " + paths(m) + " is less than or equal to zero at n=" + (n-1) + " (t=" + ((n-1)*dt) + ")");
            end
        end
    end
    paths = paths(Xn_1 > 0);
end

