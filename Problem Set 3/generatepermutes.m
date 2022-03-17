function [permutes] = generatepermutes(n)
% Generates all random walks of H or T coin flips (up or down stock
% movement) for a given number of time periods
    permutes = dec2bin(0:2^n-1) - '0';
end