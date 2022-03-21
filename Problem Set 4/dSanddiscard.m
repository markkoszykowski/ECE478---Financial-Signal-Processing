function [S, dW, paths] = dSanddiscard(S, n, beta, gamma, dt, dW, paths)
    dS = beta .* dt + gamma .* dW(:, n-1);
    Sn_1 = dS + S(:, n-1);

    S(:, n) = Sn_1 .* (Sn_1 > 0);

    S = S((Sn_1 > 0), :);
    dW = dW((Sn_1 > 0), :);
    if any(Sn_1 <= 0, "all")
        for m = find(Sn_1 <= 0).'
            disp("Path " + paths(m) + " is less than or equal to zero at n=" + (i-1) + " (t=" + ((i-1)*delta) + ")");
        end
    end
    paths = paths(Sn_1 > 0);
end

