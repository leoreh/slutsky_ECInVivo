figure

idx1 = [8 13];
idx2 = [12 13];

% manual scoring
sb1 = subplot(3, 2, 1);
plot((1 : length(sig)) / fs / 60, sig)
hold on
axis tight
Y = ylim;
fill([bs.manStamps fliplr(bs.manStamps)], [Y(1) Y(1) Y(2) Y(2)],...
    'k', 'FaceAlpha', 0.4,  'EdgeAlpha', 0);
ylabel('Voltage [mV]')
xlim([idx1])
set(gca, 'TickLength', [0 0])
box off
title('Manual Scoring')

sb1 = subplot(3, 2, 2);
plot((1 : length(sig)) / fs / 60, sig)
hold on
axis tight
Y = ylim;
fill([bs.manStamps fliplr(bs.manStamps)], [Y(1) Y(1) Y(2) Y(2)],...
    'k', 'FaceAlpha', 0.4,  'EdgeAlpha', 0);
ylabel('Voltage [mV]')
xlim(idx2)
set(gca, 'TickLength', [0 0])
box off
title('Manual Scoring')

% 1 s bins
sb1 = subplot(3, 2, 5);
plot((1 : length(sig)) / fs / 60, sig)
hold on
axis tight
Y = ylim;
fill([x.stamps fliplr(x.stamps)] / fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
    'k', 'FaceAlpha', 0.4,  'EdgeAlpha', 0);
ylabel('Voltage [mV]')
xlim(idx1)
set(gca, 'TickLength', [0 0])
box off
title('1 s minimum')

sb1 = subplot(3, 2, 6);
plot((1 : length(sig)) / fs / 60, sig)
hold on
axis tight
Y = ylim;
fill([x.stamps fliplr(x.stamps)] / fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
    'k', 'FaceAlpha', 0.4,  'EdgeAlpha', 0);
ylabel('Voltage [mV]')
xlim(idx2)
set(gca, 'TickLength', [0 0])
box off
title('1 s minimum')


% 2 s bins
sb1 = subplot(3, 2, 3);
plot((1 : length(sig)) / fs / 60, sig)
hold on
axis tight
Y = ylim;
fill([bs.stamps fliplr(bs.stamps)] / fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
    'k', 'FaceAlpha', 0.4,  'EdgeAlpha', 0);
ylabel('Voltage [mV]')
xlim(idx1)
set(gca, 'TickLength', [0 0])
box off
title('2 s minimum')

sb1 = subplot(3, 2, 4);
plot((1 : length(sig)) / fs / 60, sig)
hold on
axis tight
Y = ylim;
fill([bs.stamps fliplr(bs.stamps)] / fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
    'k', 'FaceAlpha', 0.4,  'EdgeAlpha', 0);
ylabel('Voltage [mV]')
xlim(idx2)
set(gca, 'TickLength', [0 0])
box off
title('2 s minimum')
suptitle(basename)


figure
subplot(1, 2, 1)
dur = (x.manStamps(:, 2) - x.manStamps(:, 1)) * 60;
histogram(dur, 15)
hold on
plot([2 2], ylim, '--k', 'LineWidth', 2)
title('short_Wt2')
xlabel('Burst duration [s]')
ylabel('Counts')

subplot(1, 2, 2)
dur = (bs.manStamps(:, 2) - bs.manStamps(:, 1)) * 60;
histogram(dur, 15)
hold on
plot([2 2], ylim, '--k', 'LineWidth', 2)
title('short_Wt8')
xlabel('Burst duration [s]')
ylabel('Counts')


if saveFig
    figname = fullfile(basepath, 'Duration');
    % print(fh, figname, '-dpdf', '-bestfit', '-painters');
    export_fig(figname, '-tif', '-transparent', '-r300')
end