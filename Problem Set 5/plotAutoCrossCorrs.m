function plotAutoCrossCorrs(rts, M, vars, ys, titles)
    figure;
    for i = 1:size(rts, 3)
        [rho1, lags1] = autocorr(rts(1, :, i), NumLags=M);
        [rho2, lags2] = autocorr(rts(2, :, i), NumLags=M);
        [rho12, lags12] = crosscorr(rts(1, :, i), rts(2, :, i), NumLags=M);
    
        subplot(size(rts, 3), 3, (i - 1)*3 + 1);
        plot(lags1, rho1);
        if exist("titles", "var")
            title(ys(1) + " of " + vars + " with {\itv_{t}}~" + titles(i));
        else
            title(ys(1) + " of " + vars);
        end
        xlabel("{\itm}");
        ylabel(ys(1));
    
        subplot(size(rts, 3), 3, (i - 1)*3 + 2);
        plot(lags2, rho2);
        if exist("titles", "var")
            title(ys(2) + " of " + vars + " with {\itv_{t}}~" + titles(i));
        else
            title(ys(2) + " of " + vars);
        end
        xlabel("{\itm}");
        ylabel(ys(2));
    
        subplot(size(rts, 3), 3, (i - 1)*3 + 3);
        plot(lags12, rho12);
        if exist("titles", "var")
            title(ys(3) + " of " + vars + " with {\itv_{t}}~" + titles(i));
        else
            title(ys(3) + " of " + vars);
        end
        xlabel("{\itm}");
        ylabel(ys(3));
    end
end

