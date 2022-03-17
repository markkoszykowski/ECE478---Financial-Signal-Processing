function [St, HT] = generatemrandompaths(S0, u, d, p, q, n, m, tolerance)
    assert(abs(1 - (p + q)) <= tolerance, "'p' and 'q' values do not reflect the probabilities of a valid Bernoulli trial");
    
    HT = rand(m, n) < p;
    St = zeros([size(HT, 1) n+1]);
    St(:, 1) = S0;
    for i = 2:size(St, 2)
        St(:, i) = u * St(:, i-1) .* HT(:, i-1) + d * St(:, i-1) .* ~HT(:, i-1);
    end
end

