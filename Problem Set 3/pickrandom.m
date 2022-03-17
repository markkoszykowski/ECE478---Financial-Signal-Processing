function [SN] = pickrandom(Sn, Pn, tolerance)
% Picks a random variable with a provided pmf
    assert(abs(sum(Pn) - 1) <= tolerance, "Probabilities are not valid.");
    assert(all(size(Sn) == size(Pn), "all"), "Values and probabilities vectors do not align");
    % generate cmf
    Cn = [0 cumsum(Pn)];
    SN = Sn(find(rand > Cn, 1, "last"));
end