function [Sn, Pn] = distribution(S0, u, d, p, q, n, tolerance)
% Returns stock values and their probabilities for a given time by
% evaluating all possible paths
    assert(abs(1 - (p + q)) <= tolerance, "'p' and 'q' values do not reflect the probabilities of a valid Bernoulli trial");

    % sort for numerically stable discrete values
    permutes = sort(generatepermutes(n), 2, "descend");

    % convert binary values to stock changes and probabilities
    ud = permutes * u + ~permutes * d;
    pq = permutes * p + ~permutes * q;
    
    dist = [S0*prod(ud, 2).' ; prod(pq, 2).'];

    % reduce vectors to unique stock values and respective probabilities
    [Sn, ~, Pn] = unique(dist(1, :));
    Pn = cell2mat(accumarray(Pn, 1:numel(Pn), [], @(x) {sum(dist(2:end, x), 2)})).';
end