function [SN] = pickrandom(Sn, Pn, tolerance)
    assert(abs(sum(Pn) - 1) <= tolerance, "Probabilities are not valid.");
    assert(all(size(Sn) == size(Pn)), "Values and probabilities vectors do not align");
    Cn = [0 cumsum(Pn)];
    SN = Sn(find(rand > Cn, 1, "last"));
end