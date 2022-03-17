function [Delta_n, Xn, Mn] = replicateonestep(V, Sn, u, d, r, K, tolerance, Xn_1H, Xn_1T)
% Computes value of asset, number of shares of underlying stock, and 
% number of shares in the money market at a time period using replicating 
% portfolio and wealth equation
    r = 1 + r;

    Sn_1H = u * Sn;
    Sn_1T = d * Sn;
    % if not provided with values from previous iteration (next time step),
    % generate them using provided value of asset function
    if ~exist("Xn_1H", "var")
        Xn_1H = V(Sn_1H, K);
    end
    if ~exist("Xn_1T", "var")
        Xn_1T = V(Sn_1T, K);
    end

    % solve system of 2 linear equations for each value going up or down
    a = Sn_1H - r * Sn;
    b = Sn_1T - r * Sn;

    Delta_n = (Xn_1H - Xn_1T) ./ (a - b);

    assert(all(abs((Xn_1H - Delta_n .* a) - (Xn_1T - Delta_n .* b)) <= tolerance, "all"), ...
        "Linear system of equations is not solvable");
    
    V(Sn, K);
    Xn = (Xn_1H - Delta_n .* a) / r;

    if any(Delta_n < 0, "all")
        disp("    Shorting stock");
    end

    Mn = Xn - Delta_n .* Sn;
    if any(Mn < 0, "all")
        disp("    Shorting money market");
    end

    % alternative method with matrices that can not be vectorized
    
    % A = [Sn_1H - r*Sn r ; Sn_1T - r*Sn r];
    % b = [Xn_1H ; Xn_1T];
    % x = A \ b;
    % Delta_n = x(1);
    % Xn = x(2);
end