function [X, dW, paths] = dXanddiscard(X, n, beta, gamma, dt, dW, paths)
    dX = beta .* dt + gamma .* dW(:, n-1);
    Xn_1 = dX + X(:, n-1);

    X(:, n) = Xn_1 .* (Xn_1 > 0);

    X = X((Xn_1 > 0), :);
    dW = dW((Xn_1 > 0), :);
    if any(Xn_1 <= 0, "all")
        for m = find(Xn_1 <= 0).'
            disp("Path " + paths(m) + " is less than or equal to zero at n=" + (n-1) + " (t=" + ((n-1)*dt) + ")");
        end
    end
    paths = paths(Xn_1 > 0);
end

