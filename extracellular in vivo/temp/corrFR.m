function corrFR(x1, y1, x2, y2, fh)

[rho1, pval1] = corr(log10(x1), y1);
[rho2, pval2] = corr(log10(x2), y2);
pyr_lf = polyfit(log10(x1), y1, 1);
y12 = polyval(pyr_lf, log10(x1));
int_lf = polyfit(log10(x2), y2, 1);
y22 = polyval(int_lf, log10(x2));

% plot
xLimit = [5 * 1e-2 5 * 1e1];

subplot(fh)
scatter(x1, y1, 20, 'b', 'filled')
hold on
scatter(x2, y2, 20, 'r', 'filled')
plot(x1, y12, '--b')
plot(x2, y22, '--r')
set(gca, 'XScale', 'log')
hold on
xlim(xLimit)
xticks([1e-2, 1e-1, 1e1])
xlabel('Mean Firing Rate [Hz]')
subtitle(sprintf('pyr: r=%.2f, p=%.2f \nint: r=%.2f, p=%.2f',...
    rho1, pval1, rho2, pval2))

end