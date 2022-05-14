function [a_LD, a_LS] = generateAR(rt, M, plotTitle)
% Returns AR FIR filter coefficients using Levinson-Durbin recursion and a
% Least-Squares approximation - same code from main code for Question 2
    [rho, lags] = autocorr(rt, NumLags=M);

    figure;
    plot(lags, rho);
    title(plotTitle);
    xlabel("{\it m}");
    ylabel("\rho({\itm})");

    C = toeplitz(rho);
    eigenvalues = eig(C);
    if any(eigenvalues < 0)
        disp(length(eigenvalues(eigenvalues < 0)) + " negative eigenvalues");
    end

    [L, D, perm] = ldl(C, "vector");

    F = zeros([M + 1 M + 1]);
    P = zeros([M + 1 1]);
    for m = 0:M
        [FPEF_a, P(m + 1), kappa] = levinson(rho, m);
        F(m + 1, :) = padarray(fliplr(FPEF_a), [0 M - m], "post");
    end
    disp("Reflection coefficients:");
    disp(kappa.');

    Pm = F * C * F.';

    assert(all(abs(diag(Pm) - P) <= eps*1e1), "FCF' diagonal is not equal to prediction error powers");

    if any(abs(F - L^-1) > eps*1e1, "all")
        disp("F matrix from Levinson-Durbin recursion is not equal to L^-1 from LDL decomposition");
    end
    if any(abs(diag(D) - P) > eps*1e1, "all")
        disp("P prediction error powers from Levinson-Durbin recursion is not equal to diagonal of D from LDL decomposition");
    end

    a_LD = fliplr(F(M + 1, :));
    a_LS = ar(rt, M, "ls").A;


    disp(newline + "AR coefficients from Levinson-Durbin recursion (M=" + M + "):");
    disp(a_LD);

    disp("AR coefficients from Least-Squares fit (M=" + M + "):");
    disp(a_LS);
    disp(newline);
end