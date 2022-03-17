function [Sn, Pn] = binomialdistribution(S0, u, d, p, q, n, tolerance)
    assert(abs(1 - (p + q)) <= tolerance, "'p' and 'q' values do not reflect the probabilities of a valid Bernoulli trial");

    x = 0:n;
    Pn = binopdf(x, n, p);
    Sn = S0 * u .^ x .* d .^ (n - x);
end