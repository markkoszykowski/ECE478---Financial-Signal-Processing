function [permutes] = generatepermutes(n)
    permutes = dec2bin(0:2^n-1) - '0';
end