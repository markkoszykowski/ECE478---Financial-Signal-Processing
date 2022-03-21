%% Problem Set 4 ECE478 Mark Koszykowski

clc;
clear;
close all;
%% Known Parameters

N = 250;
delta = 0.01;

t = delta * (0:N);
T = t(end);

V = @(S, K) (S > K) .* (S - K);
D = @(r, t) 1 / (1 + r) ^ t;

tolerance = 1;
tolerance = eps * 10 ^ tolerance;
%% Check Weiner Process

W = weiner(N, delta, 1000);

figure;
for i = 1:size(W, 1)
    plot(t, W(i, :));
    hold on;
end
hold off;
title("Weiner Process (\it{\mu}=" + mean(W(:, end)) + ", \it{\sigma}^2=" + var(W(:, end)) + ")");
xlabel("\it{t}");
ylabel("\it{W}(\it{t})");
%% 1

S0 = 1;
alpha = 0.1;
sigma = 0.2;
r = 0.05;

% a
theta = (alpha  - r) / sigma;

M = 1000;

dt = delta;
assert(all(abs(diff(t) - dt) <= tolerance, "all"), "'dt' does not match 'delta'")

W = weiner(N, delta, M);
dW = diff(W, 1, 2);

W_tilde = weiner(N, delta, M);
dW_tilde = diff(W_tilde, 1, 2);

% b
S = zeros(M, N+1);
S_rn = zeros(M, N+1);
S(:, 1) = S0;
S_rn(:, 1) = S0;

S_remaining = 1:M;
S_rn_remaining = 1:M;
for n = 1:N
    [S, dW, S_remaining] = dXanddiscard(S, n+1, alpha*S(:, n), sigma*S(:, n), ...
        dt, dW, S_remaining);

    [S_rn, dW_tilde, S_rn_remaining] = dXanddiscard(S_rn, n+1, r*S_rn(:, n), ...
        sigma*S_rn(:, n), dt, dW_tilde, S_rn_remaining);
end

figure;
subplot(1, 2, 1);
for i = 1:size(S, 1)
    plot(t, S(i, :));
    hold on;
end
hold off;
xlim([min(t) max(t)]);
title("\it{S}_{\it{t}} (\it{S}(0)=" + S0 + ", \it{\alpha}=" + alpha + ...
    ", \it{\sigma}=" + sigma + ")");
xlabel("\it{t}");
ylabel("\it{S}(\it{t})");

subplot(1, 2, 2);
for i = 1:size(S_rn, 1)
    plot(t, S_rn(i, :));
    hold on;
end
hold off;
xlim([min(t) max(t)]);
title("\it{S}_{\it{t}} Risk Neutral (\Theta=" + theta + ", \it{r}=" + r + ")");
xlabel("\it{t}");
ylabel("\it{S}(\it{t})");

% c
SN = mean(S(:, N + 1));
SN_2 = mean(S(:, N/2 + 1));
disp(table(M, SN, SN_2));

% d
K = SN;

SN_2 = 0:0.01:(1.25*K);
VN_2 = BSM(T - t(N/2 + 1), SN_2, K, r, sigma);

figure;
plot(SN_2, VN_2);
xlim([min(SN_2) max(SN_2)]);
title("Value of European Call Option (\it{K}=" + K + ")");
xlabel("\it{S}[\it{N}/2]");
ylabel("\it{V}[\it{N}/2]");

% e
i = (1:10).';

Si_N_2 = S(i, N/2 + 1);
Vi_N_2_BSM = BSM(T - t(N/2 + 1), Si_N_2, K, r, sigma);

W_tilde = weiner(N - N/2, delta, length(Si_N_2)*M);
dW_tilde = diff(W_tilde, 1, 2);

S_N_2_N = zeros(length(Si_N_2)*M, N/2 + 1);
for j = i.'
    S_N_2_N(((j-1)*1000 + 1):(j*1000), 1) = Si_N_2(j);
end

S_N_2_N_remaining = 1:(length(i)*M);
for n = 1:(N - N/2)
    [S_N_2_N, dW_tilde, S_N_2_N_remaining] = dXanddiscard(S_N_2_N, n+1, ...
        r*S_N_2_N(:, n), sigma*S_N_2_N(:, n), dt, dW_tilde, S_N_2_N_remaining);
end

Vi_N_2_MC = zeros(size(Si_N_2));
for j = i.'
    Vi_N_2_MC(j) = D(r, T - t(N/2 + 1)) * mean(V(S_N_2_N(((j-1)*1000 + 1):(j*1000), end), K));
end

figure;
plot(Si_N_2, Vi_N_2_BSM, "r*", Si_N_2, Vi_N_2_MC, "bo");
legend(["BSM" "Monte Carlo"], "Location", "northwest");
title("Value of V at \it{N/2} of European Call Option");
xlabel("\it{S}[\it{N}/2]");
ylabel("\it{V}[\it{N}/2]");

disp(table(i, Si_N_2, Vi_N_2_BSM, Vi_N_2_MC));
%% 2

E_Rt = @(alpha, beta, r, t) (exp(-beta * t) * r + (alpha / beta) * (1 - exp(-beta * t)));
var_Rt = @(alpha, beta, sigma, r, t) ((sigma ^ 2 / beta) * r * (exp(-beta * t) - exp(-2 * beta * t)) + ...
    ((alpha * sigma ^ 2) / (2 * beta ^ 2)) * (1 - 2 * exp(-beta * t) - exp(-2 * beta * t)));


% a
beta = 1;
alpha = 0.10 * beta;
r = 0.05;
sigma = 0.5;

assert(r > 0, "Initial condition was not met for Cox-Ingersoll-Ross Interest Rate Model");

M = 1000;

expected_successes = 50/1000;

total_M = M / expected_successes;

delta = 0.01;

T = 10;
t = 0:delta:T;

dt = delta;
assert(all(abs(diff(t) - dt) <= tolerance, "all"), "'dt' does not match 'delta'")

N = T / delta;

W = weiner(N, delta, total_M);
dW = diff(W, 1, 2);

R = zeros(total_M, N+1);
R(:, 1) = r;

R_remaining = 1:total_M;
for n = 1:N
    [R, dW, R_remaining] = dXanddiscard(R, n+1, alpha - beta*R(:, n), ...
        sigma*sqrt(R(:, n)), dt, dW, R_remaining);
end

figure;
subplot(1, 2, 1);
for i = 1:size(R, 1)
    plot(t, R(i, :));
    hold on;
end
hold off;
xlim([min(t) max(t)]);
title("Cox-Ingersoll-Ross Interest Rate Model (\it{R}(0)=" + r + ...
    ", \it{\alpha}=" + alpha + ", \it{\beta}=" + beta + ", \it{\sigma}=" + sigma + ")");
xlabel("\it{t}");
ylabel("\it{R}(\it{t})");

disp(newline + "Cox-Ingersoll-Ross Interest Rate Model discretized success rate: " + ...
    length(R_remaining) + "/" + total_M + " (" + length(R_remaining)/total_M + ")");

% b
i = 1:10;

subplot(1, 2, 2);
for j = i
    plot(t, R(j, :));
    hold on;
end
hold off;
xlim([min(t) max(t)]);
title("Cox-Ingersoll-Ross Interest Rate Model using first " + max(i) + ...
    " paths (\it{R}(0)=" + r + ", \it{\alpha}=" + alpha + ", \it{\beta}=" + ...
    beta + ", \it{\sigma}=" + sigma + ")");
xlabel("\it{t}");
ylabel("\it{R}(\it{t})");

% c
ts = [1 ; 10];

E_Rt_exp = zeros(size(ts));
var_Rt_exp = zeros(size(ts));
E_Rt_MC = zeros(size(ts));
var_Rt_MC = zeros(size(ts));
for tau = ts.'
    E_Rt_exp(ts == tau) = E_Rt(alpha, beta, r, tau);
    var_Rt_exp(ts == tau) = var_Rt(alpha, beta, sigma, r, tau);

    E_Rt_MC(ts == tau) = mean(R(:, t == tau));
    var_Rt_MC(ts == tau) = var(R(:, t == tau));
end

disp(table(ts, E_Rt_exp, E_Rt_MC, var_Rt_exp, var_Rt_MC));

figure;
subplot(1, 2, 1);
for i = 1:size(W, 1)
    plot(t, W(i, :));
    hold on;
end
hold off;
title("All Weiner Processes Used(\it{\mu}=" + mean(W(:, end)) + ", \it{\sigma}^2=" + var(W(:, end)) + ")");
xlabel("\it{t}");
ylabel("\it{W}(\it{t})");

subplot(1, 2, 2);
for i = R_remaining
    plot(t, W(i, :));
    hold on;
end
hold off;
title("Surviving Weiner Processes (\it{\mu}=" + mean(W(R_remaining, end)) + ", \it{\sigma}^2=" + var(W(R_remaining, end)) + ")");
xlabel("\it{t}");
ylabel("\it{W}(\it{t})");