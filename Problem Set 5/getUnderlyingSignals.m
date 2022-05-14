function [rts, at, bt, ut, vt, sigma0] = getUnderlyingSignals(delta, G, A, b, a, nu, N, plotSignals)
% Returns set of equations and underlying random signals that generate them
% for equal given in Problem Set
    % first generate u(t)
    ut = filter(b, a, randn([1 N]));

    % construct a simple filter to generate a(t) from u(t)
    ba = [1];
    aa = [1 -G(1, 1)];

    at = filter(ba, aa, delta(1) + ut);

    sigma0 = sqrt((1 / 16) * mean(at.^2));

    % generate random signal v(t) dependent on a(t)
    vt = [normrnd(0, sigma0, [1 N]); trnd(nu, [1 N]) / (sqrt(nu / (nu - 2)) / sigma0)];

    % construct simple filter to generate b(t) from a(t) and v(t)
    bb = [1];
    ab = [1 -G(2, 2)];

    bt = zeros([2 N]);
    for i = 1:size(bt, 1)
        bt(i, :) = filter(bb, ab, delta(2) + G(2, 1) * [0 at(1:end - 1)] + vt(i, :));
    end

    % store combinations of r1 and r2 in matrix for later usage
    rts = zeros([2 N size(bt, 1)]);
    for i = 1:size(bt, 1)
        rts(:, :, i) = A * [at; bt(i, :)];
    end

    if exist("plotSignals", "var") && plotSignals == 1
        figure;
        subplot(2, 4, 1);
        plot(ut);
        title("{\itu_{t}}");
        xlabel("{\itt}");
        ylabel("{\itu_{t}}");

        subplot(2, 4, 2);
        plot(vt(1, :));
        title("{\itv_{t}}~N(0,\sigma_{0}^{2}) (\sigma_{0}=" + sigma0 + ")");
        xlabel("{\itt}");
        ylabel("{\itv_{t}}");

        subplot(2, 4, 3);
        plot(at, "DisplayName", "{\ita_{t}}");
        hold on;
        plot(bt(1, :), "DisplayName", "{\itb_{t}}");
        title("{\ita_{t}} and {\itb_{t}} with {\itv_{t}}~N(0,\sigma_{0}^{2}) (\sigma_{0}=" + sigma0 + ")");
        legend("Location", "Best");
        xlabel("{\itt}");
        ylabel("{\itv_{t}}");

        subplot(2, 4, 4);
        plot(rts(1, :, 1), "DisplayName", "{\itr_{1t}}");
        hold on;
        plot(rts(2, :, 1), "DisplayName", "{\itr_{2t}}");
        title("{\itr_{1t}} and {\itr_{2t}} with {\itv_{t}}~N(0,\sigma_{0}^{2}) (\sigma_{0}=" + sigma0 + ")");
        legend("Location", "Best");
        xlabel("{\itt}");
        ylabel("{\itr_{t}}");

        subplot(2, 4, 5);
        plot(ut);
        title("{\itu_{t}}");
        xlabel("{\itt}");
        ylabel("{\itu_{t}}");

        subplot(2, 4, 6);
        plot(vt(2, :));
        title("{\itv_{t}}~Students' T (\nu=" + nu + ")");
        xlabel("{\itt}");
        ylabel("{\itv_{t}}");

        subplot(2, 4, 7);
        plot(at, "DisplayName", "{\ita_{t}}");
        hold on;
        plot(bt(2, :), "DisplayName", "{\itb_{t}}");
        title("{\ita_{t}} and {\itb_{t}} with {\itv_{t}}~Students' T (\nu=" + nu + ")");
        legend("Location", "Best");
        xlabel("{\itt}");
        ylabel("{\itv_{t}}");

        subplot(2, 4, 8);
        plot(rts(1, :, 2), "DisplayName", "{\itr_{1t}}");
        hold on;
        plot(rts(2, :, 2), "DisplayName", "{\itr_{2t}}");
        title("{\itr_{1t}} and {\itr_{2t}} with {\itv_{t}}~Students' T (\nu=" + nu + ")");
        legend("Location", "Best");
        xlabel("{\itt}");
        ylabel("{\itr_{t}}");
    end
end