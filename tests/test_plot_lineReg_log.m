function test_plot_lineReg_log()
    % TEST_PLOT_LINEREG_LOG Verifies plot_lineReg on linear/log scales.

    hFig = figure('Name', 'plot_lineReg Test', 'Position', [100 100 1200 800]);

    %% 1. Linear-Linear
    % Expected: y = 2x + 5
    subplot(2, 2, 1);
    x = 1:10;
    y = 2*x + 5 + randn(size(x))*0.1;
    plot(x, y, 'ko'); 
    title('Linear-Linear (y=2x+5)');
    [h, stats] = plot_lineReg(x, y);
    xlabel(sprintf('Slope=%.2f (Exp: 2.00)', stats.slope));

    %% 2. Log-Log (Power Law)
    % Expected: y = 100 * x^-1.5  => log(y) = -1.5*log(x) + 2
    subplot(2, 2, 2);
    x = logspace(0, 2, 20);
    y = 100 .* x.^(-1.5) .* 10.^(randn(size(x))*0.05);
    loglog(x, y, 'ko');
    title('Log-Log (Pow: y=100 x^{-1.5})');
    [h, stats] = plot_lineReg(x, y);
    xlabel(sprintf('Slope=%.2f (Exp: -1.50)', stats.slope));
    ylim([0.1, 200]); xlim([0.8, 120]);

    %% 3. Semilog-Y (Exponential)
    % Expected: y = 10 * 10^(0.5x) => log10(y) = 0.5x + 1
    subplot(2, 2, 3);
    x = 0:0.5:5;
    y = 10 .* 10.^(0.5*x) .* 10.^(randn(size(x))*0.05);
    plot(x, y, 'ko');
    set(gca, 'YScale', 'log');
    title('Semilog-Y (Exp: y=10^{0.5x+1})');
    [h, stats] = plot_lineReg(x, y);
    xlabel(sprintf('Slope=%.2f (Exp: 0.50)', stats.slope));

    %% 4. Semilog-X (Logarithmic)
    % Expected: y = 2*log10(x) + 3
    subplot(2, 2, 4);
    x = logspace(0, 2, 20);
    y = 2*log10(x) + 3 + randn(size(x))*0.1;
    plot(x, y, 'ko');
    set(gca, 'XScale', 'log');
    title('Semilog-X (Log: y=2log(x)+3)');
    [h, stats] = plot_lineReg(x, y);
    xlabel(sprintf('Slope=%.2f (Exp: 2.00)', stats.slope));

    disp('Test completed. Check figure for straight fit lines and correct slopes.');
end
