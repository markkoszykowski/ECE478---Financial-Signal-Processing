function [St, HT] = generatepaths(S0, u, d, n)
% Generates all random walks of H or T coin flips (up or down stock
% movement) for a given number of time periods and stock values over walks
    % get all possible random walks of length n
    HT = generatepermutes(n);
    % compute stock values over time intervals
    St = zeros([size(HT, 1) n+1]);
    St(:, 1) = S0;
    for i = 2:size(St, 2)
        St(:, i) = u * St(:, i-1) .* HT(:, i-1) + d * St(:, i-1) .* ~HT(:, i-1);
    end
end