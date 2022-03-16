%% Problem Set 3 ECE478 Mark Koszykowski

clc;
clear;
close all;

tolerance = 1;

tolerance = eps * 10 ^ tolerance;
%% Known Parameters

S0 = 1;
u = 1.07;
d = 1.02;
r = 0.03;

V = @(S, K) (S > K) .* (S - K);
tilde = @(X, r, n) 1 / (1 + r) ^ n * X;
%% 1


disp("1)");
assert((d < 1+r) & (1+r < u), "Parameters do not satisfy no-artbitrage condition");

p = 2/3;
q = 1-p;

n = 10;

[Sn, Pn] = distribution(S0, u, d, p, q, n, tolerance);

disp(table(Sn.', Pn.', 'VariableNames', {'Sn', 'Pn'}));

figure;
plot(Sn, Pn);
xticks(Sn);
xtickangle(45);
title("Distribution of \it{S_{n}}");
xlabel("\it{S_{n}}");
ylabel("\it{P}(\it{S_{n}})");

Rn = log(Sn / S0);

Yn = 0:n;
Pyn = binopdf(0:n, n, p);

coeffs = polyfit(Yn, Rn, 1);
cn = coeffs(1);
dn = coeffs(2);

assert(all(abs(Pn - Pyn) <= tolerance & ...
    abs(Rn - (cn*Yn + dn)) <= tolerance), "Distribution of St is not log-binomial");


%% 2


% a
disp(newline);disp("2a)");
K = 1.42;
p_tilde = ((1 + r) - d) / (u - d);
q_tilde = (u - (1 + r)) / (u - d);

[~, Pn_tilde] = distribution(S0, u, d, p_tilde, q_tilde, n, tolerance);
disp(table(Sn.', Pn_tilde.', 'VariableNames', {'Sn', 'Pn_tilde'}));

V_tilde = tilde(V(Sn, K), r, n);
V0 = sum(V_tilde .* Pn_tilde);
disp(table(S0, u, d, r, p, q, p_tilde, q_tilde, n, K, V0));

% b
disp(newline);disp("2b)");
Sn = pickrandom(Sn, Pn, tolerance);

[Delta_n, Vn, Mn] = replicateonestep(V, Sn, u, d, r, K, tolerance);
assert(abs(Vn - (sum([V(u * Sn, K) V(d * Sn, K)] .* [p_tilde q_tilde]) / (1+r))) <= tolerance, ...
    "Wealth equation did not work");

disp(table(Sn, Delta_n, Vn, Mn))

% c
disp(newline);disp("2c)");
u = 1.10;
d = 1.01;
r = 0.05;

n = 5;

p_tilde = ((1 + r) - d) / (u - d);
q_tilde = (u - (1 + r)) / (u - d);

[SN, PN_tilde] = distribution(S0, u, d, p_tilde, q_tilde, n, tolerance);

K = sum(SN .* PN_tilde);

[Sn, omega_n] = generatepaths(S0, u, d, n);

Delta_n_omega = zeros(size(omega_n));
Mn = zeros(size(omega_n));
Vn_omega = zeros(size(Sn));

Vn_omega(:, end) = V(Sn(:, end), K);
for i = size(Delta_n_omega, 2):-1:1
    fprintf("n=%d:\n", i-1);
    skip = 2^(size(Vn_omega, 2) - i - 1);
    [Delta_n_omega(:, i), Vn_omega(:, i), Mn(:, i)] = replicateonestep(V, Sn(:, i), u, d, r, K, tolerance, ...
        repelem(Vn_omega(1+skip:skip*2:end, i+1), skip*2, 1), repelem(Vn_omega(1:skip*2:end, i+1), skip*2, 1));
end

assert(all(abs(Vn_omega(:, 1) - (sum(V(SN, K) .* PN_tilde) / (1+r)^n)) <= tolerance, "all"), ...
    "Wealth equation did not work");

V0 = Vn_omega(:, 1);
omega_n = char('H' * omega_n + 'T' * ~omega_n);

disp(table(omega_n, Delta_n_omega, Mn, V0));