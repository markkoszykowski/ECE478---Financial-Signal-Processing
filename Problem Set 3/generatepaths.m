function [St, HT] = generatepaths(S0, u, d, n)
    HT = generatepermutes(n);
    St = zeros([size(HT, 1) n+1]);
    St(:, 1) = S0;
    for i = 2:size(St, 2)
        St(:, i) = u * St(:, i-1) .* HT(:, i-1) + d * St(:, i-1) .* ~HT(:, i-1);
    end
end