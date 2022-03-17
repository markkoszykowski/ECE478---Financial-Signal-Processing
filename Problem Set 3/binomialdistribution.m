function [Sn, Pn] = binomialdistribution(S0, u, d, p, q, n, tolerance)
% Returns stock values and their probabilities for a given time using a
% binomial distribution
    assert(abs(1 - (p + q)) <= tolerance, "'p' and 'q' values do not reflect the probabilities of a valid Bernoulli trial");

    x = 0:n;
    Pn = binopdf(x, n, p);
    % compute stock values from binomial random values
    Sn = S0 * u .^ x .* d .^ (n - x);
end