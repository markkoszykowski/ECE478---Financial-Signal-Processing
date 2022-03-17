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

N = 10;

[Sn, Pn] = distribution(S0, u, d, p, q, N, tolerance);

disp(table(Sn.', Pn.', 'VariableNames', {'Sn', 'Pn'}));

figure;
plot(Sn, Pn);
xticks(Sn);
xtickangle(45);
title("Distribution of \it{S_{n}}");
xlabel("\it{S_{n}}");
ylabel("\it{P}(\it{S_{n}})");

Rn = log(Sn / S0);

Yn = 0:N;
Pyn = binopdf(0:N, N, p);

coeffs = polyfit(Yn, Rn, 1);
cn = coeffs(1);
dn = coeffs(2);

assert(all(abs(Pn - Pyn) <= tolerance &abs(Rn - (cn*Yn + dn)) <= tolerance, "all"), ...
    "Distribution of St is not log-binomial");


%% 2


% a
disp(newline);disp("2a)");
K = 1.42;
p_tilde = ((1 + r) - d) / (u - d);
q_tilde = (u - (1 + r)) / (u - d);

[~, Pn_tilde] = distribution(S0, u, d, p_tilde, q_tilde, N, tolerance);
disp(table(Sn.', Pn_tilde.', 'VariableNames', {'Sn', 'Pn_tilde'}));

V0 = sum(tilde(V(Sn, K), r, N) .* Pn_tilde);
disp(table(S0, u, d, r, p, q, p_tilde, q_tilde, N, K, V0));

% b
disp(newline);disp("2b)");
Sn = pickrandom(Sn, Pn, tolerance);

[Delta_n, Vn, Mn] = replicateonestep(V, Sn, u, d, r, K, tolerance);
assert(abs(Vn - sum(tilde([V(u*Sn, K) V(d*Sn, K)], r, 1) .* [p_tilde q_tilde])) <= tolerance, ...
    "Wealth equation did not work");

disp(table(Sn, Delta_n, Vn, Mn))

% c
disp(newline);disp("2c)");
u = 1.10;
d = 1.01;
r = 0.05;

N = 5;

p_tilde = ((1 + r) - d) / (u - d);
q_tilde = (u - (1 + r)) / (u - d);

[SN, PN_tilde] = distribution(S0, u, d, p_tilde, q_tilde, N, tolerance);

K = sum(SN .* PN_tilde);

[Sn, omega_n] = generatepaths(S0, u, d, N);

Delta_n_omega = zeros(size(omega_n));
Vn_omega = zeros(size(Sn));
Mn_omega = zeros(size(omega_n));

Vn_omega(:, end) = V(Sn(:, end), K);
for i = size(Delta_n_omega, 2):-1:1
    fprintf("n=%d:\n", i-1);
    skip = 2^(size(Vn_omega, 2) - i - 1);
    [Delta_n_omega(:, i), Vn_omega(:, i), Mn_omega(:, i)] = replicateonestep(V, Sn(:, i), u, d, r, K, tolerance, ...
        repelem(Vn_omega(1+skip:skip*2:end, i+1), skip*2, 1), repelem(Vn_omega(1:skip*2:end, i+1), skip*2, 1));
end

assert(all(abs(Vn_omega(:, 1) - sum(tilde(V(SN, K), r, N) .* PN_tilde)) <= tolerance, "all"), ...
    "Wealth equation did not work");

V0 = Vn_omega(:, 1);
omega_n = char('H' * omega_n + 'T' * ~omega_n);

disp(table(omega_n, Delta_n_omega, Mn_omega, V0));

%% 3


% a
disp(newline);disp("3a)");

V0 = sum(tilde(V(SN, K), r, N) .* PN_tilde);

M = [1 5 10 32];

estimates = table;
for m = M
    [S0_est, V0_est] = wrapper(V, S0, u, d, r, K, N, m, tolerance);
    estimates = [estimates ; table(m, S0, S0_est, u, d, r, K, N, V0_est, V0)];
end
disp(estimates);

% b
disp(newline);disp("3b)");
u = 1 + 5e-3;
d = 1 + 1e-4;
r = 1e-3;

N = 100;

p_tilde = ((1 + r) - d) / (u - d);
q_tilde = (u - (1 + r)) / (u - d);

[SN, PN_tilde] = binomialdistribution(S0, u, d, p_tilde, q_tilde, N, tolerance);

K = sum(SN .* PN_tilde);

V0 = sum(tilde(V(SN, K), r, N) .* PN_tilde);

n = 10;

[SnH, PnH_tilde] = binomialdistribution(S0*u^n, u, d, p_tilde, q_tilde, N-n, tolerance);
[SnT, PnT_tilde] = binomialdistribution(S0*d^n, u, d, p_tilde, q_tilde, N-n, tolerance);

VnH = sum(tilde(V(SnH, K), r, N-n) .* PnH_tilde);
VnT = sum(tilde(V(SnT, K), r, N-n) .* PnT_tilde);

p0 = [0.9 1.1];
S0p0s = zeros(size(p0));
V0p0s = zeros(size(p0));
for p_tilde_scalar = p0
    [SNp0, PNp0] = binomialdistribution(S0, u, d, p_tilde_scalar*p_tilde, ...
        1 - (p_tilde_scalar*p_tilde), N, tolerance);
    S0p0s(p0 == p_tilde_scalar) = sum(tilde(SNp0, r, N) .* PNp0);

    V0p0s(p0 == p_tilde_scalar) = sum(tilde(V(SNp0, K), r, N) .* PNp0);
end

M = [100 1000 10000 100000];

estimates = table;
Vn_estimates = table;
actual_p_estimates = table;
for m = M
    [S0_est, V0_est] = wrapper(V, S0, u, d, r, K, N, m, tolerance);
    estimates = [estimates ; table(m, S0, S0_est, u, d, r, K, N, V0_est, V0)];

    [~, VnH_est] = wrapper(V, S0*u^n, u, d, r, K, N-n, m, tolerance);
    [~, VnT_est] = wrapper(V, S0*d^n, u, d, r, K, N-n, m, tolerance);
    Vn_estimates = [Vn_estimates ; table(m, S0, u, d, r, K, N, n, VnH_est, VnH, VnT_est, VnT)];

    for p_tilde_scalar = p0
        p = p_tilde_scalar * p_tilde;
        S0p0 = S0p0s(p0 == p_tilde_scalar);
        V0p0 = V0p0s(p0 == p_tilde_scalar);
        [S0_est_p0, V0_est_p0] = wrapper(V, S0, u, d, r, K, N, m, tolerance, p);
        actual_p_estimates = [actual_p_estimates ; table(m, p_tilde_scalar, p_tilde, p, ...
            S0, S0_est_p0, S0p0, u, d, r, K, N, V0_est_p0, V0p0, V0)];
    end
end
disp(estimates);
disp(Vn_estimates);
disp(actual_p_estimates);