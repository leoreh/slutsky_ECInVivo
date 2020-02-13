% plotAnestatesSU

marg = 0.05;
wvstamps = linspace(-marg, marg, floor(marg * lfp.fs) * 2 + 1);
sig = double(lfp.data(:, ch));
tstamps = (1 : length(sig)) / lfp.fs / 60;

fh = figure('Visible', 'on');
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

subplot(4, 2, 1 : 2);
plot(tstamps, sig, 'k')
axis tight
hold on
plot(xlim, [iis.thr(2) iis.thr(2)], '--r')
ylabel('Voltage [uV]')
y = ylim;
ylim([y(1) 2000])
yyaxis right
plot(iis.cents / lfp.fs / 60, iis.rate, 'r', 'LineWidth', 3);
p1.Color(4) = 0.3;
ylim([0 1])
ylabel('Rate [spikes / bin]')
legend({'Raw', 'IIS thr', 'IIS rate'}, 'Location', 'northwest')
axis tight
set(gca, 'TickLength', [0 0])
box off
title('raw and IIS rate')
addLns('lns', info.lns, 'lbs', info.labels, 'ax', 'x')

subplot(4, 2, 3 : 4);
plotFRtime('fr', fr, 'units', true, 'spktime', spikes.times,...
    'avg', false, 'lns', info.lns, 'lbs', info.labels,...
    'raster', false, 'saveFig', false);

subplot(4, 2, 5 : 6)
plot(bs.cents / lfp.fs / 60, bs.bsr, 'k', 'LineWidth', 2)
hold on
plot(bs.cents / lfp.fs / 60, ep.dband, 'b', 'LineWidth', 2)
ylabel('a.u.')
legend({'BSR', 'Delta', 'IIS rate'}, 'Location', 'northeast')
axis tight
set(gca, 'TickLength', [0 0])
box off
addLns('lns', info.lns, 'lbs', info.labels, 'ax', 'x')
title('Anesthesia State')

% zoom in raw, bs, and iis
subplot(4, 2, 7);
minmarg = 1.5;
midsig = round(length(sig) / 2);
idx = round(midsig - minmarg * lfp.fs * 60 : midsig + minmarg * lfp.fs * 60);
idx2 = iis.peakPos > idx(1) & iis.peakPos < idx(end);
plot(tstamps(idx), sig(idx), 'k')
axis tight
hold on
x = xlim;
plot(x, [iis.thr(2) iis.thr(2)], '--r')
scatter(iis.peakPos(idx2) / lfp.fs / 60,...
    iis.peakPower(idx2), '*');
bsstamps = RestrictInts(bs.stamps, [idx(1) idx(end)]);
Y = ylim;
if ~isempty(bsstamps)
    fill([bsstamps fliplr(bsstamps)] / lfp.fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
        'k', 'FaceAlpha', 0.25,  'EdgeAlpha', 0);
end
ylabel('Voltage [uV]')
xlabel('Time [m]')
xticks([ceil(midsig / lfp.fs / 60 - minmarg), floor(midsig / lfp.fs / 60 + minmarg)])
set(gca, 'TickLength', [0 0])
box off
title('IIS')

% iis waveforms
subplot(4, 2, 8)
plot(wvstamps * 1000, iis.wv)
ylabel('Voltage [mV]')
xlabel('Time [ms]')
axis tight
xticks([-marg, 0, marg] * 1000);
set(gca, 'TickLength', [0 0])
box off
title('IIS waveform')

% mean + std waveform
axes('Position',[.571 .11 .15 .1])
box on
stdshade(iis.wv, 0.5, 'k', wvstamps)
axis tight
set(gca, 'TickLength', [0 0], 'YTickLabel', [], 'XTickLabel', [],...
    'XColor', 'none', 'YColor', 'none', 'Color', 'none')
title(sprintf('n = %d', size(iis.wv, 1)));
box off


[~, basename] = fileparts(basepath);
figname = [basename '_anestates'];
export_fig(figname, '-tif', '-transparent')