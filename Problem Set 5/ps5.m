%% Problem Set 5 ECE478 Mark Koszykowski

clc;
clear;
close all;
%% 1


N = 1e6;

heavyTailDistributions = containers.Map();

normalCauchyScale = @(alpha) (1 / pi) * atan(1 / alpha) + (1 / 2) ...
     - ((1 / pi) * atan(-1 / alpha) + (1 / 2)) ...
     - (normcdf(1) - normcdf(-1));
alpha = fzero(normalCauchyScale, 0.5);

nus = [5 10];
heavyTailDistributions("Normal") = randn([1 N]);
heavyTailDistributions("Cauchy (\alpha=" + alpha + ")") = alpha * trnd(1, [1 N]);
for nu = nus
    heavyTailDistributions("Students' T (\nu=" + nu + ")") = trnd(nu, [1 N]) / sqrt(nu / (nu - 2));
end

figure;
i = 1;
fractions = zeros([1 length(heavyTailDistributions.keys)]);
for distribution = heavyTailDistributions.keys
    distribution = string(distribution);
    samples = heavyTailDistributions(distribution);

    subplot(2, 2, i)
    histogram(samples);
    title(distribution + " Distribution");
    xlabel("{\it x}");
    ylabel("{\it f(x)}");

    fractions(i) = length(samples(abs(samples) > 4)) / length(samples);
    i = i + 1;
end

figure;
bar(fractions);
title("Fraction of |X| > 4 for Different Distributions (N=" + N + ")");
ylabel("Fraction");
xticklabels(heavyTailDistributions.keys);

text(1:length(fractions), fractions, num2str(fractions'), ...
    "vert", "bottom", "horiz", "center"); 
box off
%% 2


