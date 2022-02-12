%% Problem Set 1 ECE478 Mark Koszykowski

clc;
clear;
close all;
%% 1

long = @(Dt) Dt;
short = @(Dt) -Dt;

call = @(St, K) (St > K).*(St - K);
put = @(St, K) (K > St).*(K - St);

digi_call = @(St, K) (St > K);
digi_put = @(St, K) (St < K);

St = 0:.001:3;

return_ylim = @(s) [min(s)-.05*(max(s)-min(s)) max(s)+.05*(max(s)-min(s))];


% a

K = 1;

long_call = long(call(St, K));
long_put = long(put(St, K));
short_call = short(call(St, K));
short_put = short(put(St, K));

figure;

subplot(2, 2, 1);
plot(St, long_call);
title("Long Call (\it{K}=" + K + ")")
xlabel("\it{S_{T}}");
ylabel("\it{V}(\it{S_{T}})");
ylim(return_ylim(long_call));

subplot(2, 2, 2);
plot(St, long_put);
title("Long Put (\it{K}=" + K + ")")
xlabel("\it{S_{T}}");
ylabel("\it{V}(\it{S_{T}})");
ylim(return_ylim(long_put));

subplot(2, 2, 3);
plot(St, short_call);
title("Short Call (\it{K}=" + K + ")")
xlabel("\it{S_{T}}");
ylabel("\it{V}(\it{S_{T}})");
ylim(return_ylim(short_call));

subplot(2, 2, 4);
plot(St, short_put);
title("Short Put (\it{K}=" + K + ")")
xlabel("\it{S_{T}}");
ylabel("\it{V}(\it{S_{T}})");
ylim(return_ylim(short_put));


% b

K = 1;

straddle_emp = call(St, K) + put(St, K);
straddle_theo = abs(St - K);

figure;

subplot(1, 2, 1);
plot(St, straddle_emp);
title("Straddle (\it{K}=" + K + ")")
xlabel("\it{S_{T}}");
ylabel("\it{V}(\it{S_{T}})");
ylim(return_ylim(straddle_emp));

subplot(1, 2, 2);
plot(St, straddle_theo);
title("\it{V}(\it{S_{T}}) = \it{|S_{T} - K|} (\it{K}=" + K + ")")
xlabel("\it{S_{T}}");
ylabel("\it{V}(\it{S_{T}})");
ylim(return_ylim(straddle_theo));

assert(isequal(straddle_emp, straddle_theo));


% c

K1 = 0.5;
K2 = 1.25;

assert(K1 < K2);

call_put_spread = long(call(St, K1)) + short(call(St, K2));

figure;
plot(St, call_put_spread);
title("Call-Put Spread (\it{K_{1}}=" + K1 + ", \it{K_{2}}=" + K2 + ")")
xlabel("\it{S_{T}}");
ylabel("\it{V}(\it{S_{T}})");
ylim(return_ylim(call_put_spread));

% d

K1 = 0.5;
K2 = 1.25;

lambda = [1/3, 1/2, 2/3];

assert(K1 < K2);
assert(all(0 < lambda) && all(lambda < 1));

K_star = lambda*K1 + (1-lambda)*K2;

figure;
for i = 1:numel(lambda)

    butterfly = long(lambda(i)*call(St, K1)) + ...
        long((1 - lambda(i))*call(St, K2)) + ...
        short(call(St, K_star(i)));

    subplot(1, 3, i);
    plot(St, butterfly);
    title("Butterfly (\it{K_{1}}=" + K1 + ...
        ", \it{K_{2}}=" + K2 + ", \it{\lambda}=" + lambda(i) + ...
        ", \it{K*}=" + K_star(i) + ")");
    xlabel("\it{S_{T}}");
    ylabel("\it{V}(\it{S_{T}})");
    ylim(return_ylim(butterfly));
end


% e

K1 = 0.5;
K2 = 1.0;
K3 = 1.5;

assert((K1 < K2) && (K2 < K3));

call_ladder = long(call(St, K1)) + short(call(St, K2)) + short(call(St, K3));

figure;
plot(St, call_ladder);
title("Call Ladder (\it{K_{1}}=" + K1 + ...
    ", \it{K_{2}}=" + K2 + ...
    ", \it{K_{3}}=" + K3 + ")")
xlabel("\it{S_{T}}");
ylabel("\it{V}(\it{S_{T}})");
ylim(return_ylim(call_ladder));


% f

K1 = 0.5;
K2 = 1.25;

assert(K1 < K2);

digi_call_spread = long(digi_call(St, K1)) + short(digi_call(St, K2));

figure;
plot(St, digi_call_spread);
title("Digital Call Spread (\it{K_{1}}=" + K1 + ...
    ", \it{K_{2}}=" + K2 + ")")
xlabel("\it{S_{T}}");
ylabel("\it{V}(\it{S_{T}})");
ylim(return_ylim(digi_call_spread));