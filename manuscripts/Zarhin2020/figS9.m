% figS9 in Zarhin et al., 2020

close all
force = false;

basepath{1} = 'F:\Daniel IIS 11.6.20\IIS\WT\New analyis';
basepath{2} = 'F:\Daniel IIS 11.6.20\IIS\APPPS1\New analysis';

% select mouse (according to file order in path)
i_wt = 7;
i_app = 27;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if force
    
    % load wt data
    cd(basepath{1})
    filename = dir('*.lfp.*');
    files = natsort({filename.name});
    nfiles = 1 : length(files);
    [~, basename] = fileparts(files{nfiles(i_wt)});
    [~, basename] = fileparts(basename);
    [bs, iis, ep] = aneStates('ch', 1, 'basepath', basepath{1},...
        'basename', basename, 'graphics', false,...
        'saveVar', false, 'saveFig', false, 'forceA', false,...
        'binsize', 30, 'smf', 7, 'thrMet', 1);
    load([fullfile(basepath{1}, basename) '.lfp.mat'])
    lfpwt = lfp;
    bswt = bs;
    iiswt = iis;
    epwt = ep;
    
    % load apps1 data
    cd(basepath{2})
    filename = dir('*.lfp.*');
    files = natsort({filename.name});
    nfiles = 1 : length(files);
    [~, basename] = fileparts(files{nfiles(i_app)});
    [~, basename] = fileparts(basename);
    [bs, iis, ep] = aneStates('ch', 1, 'basepath', basepath{2},...
        'basename', basename, 'graphics', false,...
        'saveVar', false, 'saveFig', false, 'forceA', false,...
        'binsize', 30, 'smf', 7, 'thrMet', 1);
    load([fullfile(basepath{2}, basename) '.lfp.mat'])

    
    % bs classification
    fs = lfp.fs;
    binsize = 2 ^ nextpow2(0.5 * fs);   % for pc1

    ibins = [1 : binsize : length(lfpwt.data)];
    ibins(end) = length(lfpwt.data);
    
    % divide signal to bins
    sigre = lfpwt.data(1 : end - (mod(length(lfpwt.data), binsize) + binsize));
    sigmat = reshape(sigre, binsize, (floor(length(lfpwt.data) / binsize) - 1));
    
    % last bin
    siglastbin = lfpwt.data(length(sigre) + 1 : length(lfpwt.data));
    
    stdvec = std(sigmat);
    stdvec = [stdvec, std(siglastbin)];
    maxvec = max(abs(sigmat));
    maxvec = [maxvec, max(abs(siglastbin))];
    freq = logspace(0, 2, 100);
    win = hann(binsize);
    [~, fff, ttt, pband] = spectrogram(lfpwt.data, win, 0, freq, fs, 'yaxis', 'psd');
    pband = 10 * log10(abs(pband));
    [~, pc1] = pca(pband', 'NumComponents', 1);
    
    % concatenate and lognorm
    varmat = [stdvec; maxvec];
    if size(varmat, 1) < size (varmat, 2)
        varmat = varmat';
    end
    varmat = log10(varmat);
    varmat = [varmat, pc1];
    
    options = statset('MaxIter', 500);
    gm = fitgmdist(varmat, 2,...
        'options', options, 'Start', 'plus', 'Replicates', 50);
    % cluster
    [gi, ~, ~, ~, mDist] = cluster(bs.gm, varmat);
    
    options = statset('MaxIter', 50);
    gm = fitgmdist(varmat(:, [1, 3]), 2,...
        'options', options, 'Start', 'plus', 'Replicates', 50);
    
    % group data
    load(fullfile('F:\Daniel IIS 11.6.20\IIS', 'as.mat'))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% general params
fh = figure('Visible', 'on');
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

Y = [-1 2]; % ylim for raw
YY = [-0.5 1]; % ylim for zoom
binsize = (2 ^ nextpow2(fs * 2));
smf = 15;
marg = 0.05;
minmarg = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wt mouse

% idx for total recording
idx1 = 5 * fs * 60 : 25 * fs * 60;

% idx for zoomin in samples
midsig = 23.5;
idx2 = round((midsig - minmarg) * fs * 60 : (midsig + minmarg) * fs * 60);

% raw
subplot(6, 2, 1);
plot(lfpwt.timestamps(idx1) / 60, lfpwt.data(idx1), 'k', 'LineWidth', 1)
hold on
set(gca, 'TickLength', [0 0], 'Color', 'none', 'XTickLabel', [],...
    'XColor', 'none')
yticks(Y)
box off
ylabel('Voltage [mV]')
ylim(Y)
idx3 = [idx2(1) idx2(end)] / fs / 60;
fill([idx3 fliplr(idx3)]', [Y(1) Y(1) Y(2) Y(2)],...
    'r', 'FaceAlpha', 0.2,  'EdgeAlpha', 0, 'HandleVisibility', 'off');
axis tight
xlabel('Time [m]')
title('WT')

% spectrogram
subplot(6, 2, 3);
specBand('sig', lfpwt.data(idx1), 'graphics', true, 'binsize', binsize,...
    'smf', smf, 'normband', true);
set(gca, 'TickLength', [0 0], 'XTickLabel', [],...
    'Color', 'none', 'XColor', 'none')
ylim([0 100])
title('')
colorbar('off');

% bsr
subplot(6, 2, 5);
[~, bsidx] = min(abs(bswt.cents - idx1(1)));
[~, bsidx(2)] = min(abs(bswt.cents - idx1(end)));
bsidx = bsidx(1) : bsidx(2);
plot(bswt.cents(bsidx) / fs / 60, bswt.bsr(bsidx), 'k', 'LineWidth', 1)
hold on
plot(bswt.cents(bsidx) / fs / 60, epwt.dband(bsidx), 'b', 'LineWidth', 1)
axis tight
ylim([0 1])
set(gca, 'TickLength', [0 0],...
    'Color', 'none')
box off
yticks([0 1])
legend({'BSR', 'Delta Power'}, 'Location', 'northwest')
xlim([idx1(1) idx1(end)] / fs / 60)
xlabel('Time [m]')
xticks(idx1(1) / 60 / fs : 5 : idx1(end) / 60 / fs)
xticklabels({'1', '5', '10', '15', '20'})

% zoom in
subplot(6, 2, 7);
idx5 = iiswt.peakPos > idx2(1) & iiswt.peakPos < idx2(end);
plot(lfpwt.timestamps(idx2) / 60, lfpwt.data(idx2), 'k')
axis tight
hold on
x = xlim;
plot(x, [iiswt.thr(2) iiswt.thr(2)], '--r')
scatter(iiswt.peakPos(idx5) / fs / 60,...
    iiswt.peakPower(idx5), '*');
bsstamps = RestrictInts(bswt.stamps, [idx2(1) - fs * 20 idx2(end) + 20 * fs]);
ylim(YY)
yticks(YY)
if ~isempty(bsstamps)
    fill([bsstamps fliplr(bsstamps)] / fs / 60, [YY(1) YY(1) YY(2) YY(2)],...
        'k', 'FaceAlpha', 0.25,  'EdgeAlpha', 0);
end
ylabel('Voltage [mV]')
xlabel('Time [m]')
xlim(idx3)
xticks(idx3)
xticklabels([17.5 19.5])
set(gca, 'TickLength', [0 0])
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% app mouse

% idx for total recording
idx1 = 15 * fs * 60 : 35 * fs * 60;

% idx for zoomin in samples
midsig = 33.5;
idx2 = round((midsig - minmarg) * fs * 60 : (midsig + minmarg) * fs * 60);

% raw
subplot(6, 2, 2);
plot(lfp.timestamps(idx1) / 60, lfp.data(idx1), 'k', 'LineWidth', 1)
hold on
set(gca, 'TickLength', [0 0], 'XTickLabel', [], 'XColor', 'none',...
    'Color', 'none', 'YColor', 'none')
box off
ylim(Y)
% yticks(Y)
idx3 = [idx2(1) idx2(end)] / fs / 60;
fill([idx3 fliplr(idx3)]', [Y(1) Y(1) Y(2) Y(2)],...
    'r', 'FaceAlpha', 0.2,  'EdgeAlpha', 0, 'HandleVisibility', 'off');
axis tight
xlabel('Time [m]')
% ylabel('LFP [mV]')
title('APP-PS1', 'Interpreter', 'none')

% spectrogram
subplot(6, 2, 4);
specBand('sig', lfp.data(idx1), 'graphics', true, 'binsize', binsize,...
    'smf', smf, 'normband', true);
set(gca, 'TickLength', [0 0], 'XTickLabel', [],...
    'Color', 'none', 'XColor', 'none', 'YColor', 'none')
ylim([0 100])
title('')

% bsr
subplot(6, 2, 6);
[~, bsidx] = min(abs(bs.cents - idx1(1)));
[~, bsidx(2)] = min(abs(bs.cents - idx1(end)));
bsidx = bsidx(1) : bsidx(2);
plot(bs.cents(bsidx) / fs / 60, bs.bsr(bsidx), 'k', 'LineWidth', 1)
hold on
plot(bs.cents(bsidx) / fs / 60, ep.dband(bsidx), 'b', 'LineWidth', 1)
set(gca, 'TickLength', [0 0], 'YTickLabel', [],...
    'Color', 'none', 'YColor', 'none')
box off
axis tight
ylim([0 1])
xlim([idx1(1) idx1(end)] / fs / 60)
xlabel('Time [m]')
xticks(idx1(1) / 60 / fs : 5 : idx1(end) / 60 / fs)
xticklabels({'1', '5', '10', '15', '20'})

% zoom in
subplot(6, 2, 8);
idx5 = iis.peakPos > idx2(1) & iis.peakPos < idx2(end);
plot(lfp.timestamps(idx2) / 60, lfp.data(idx2), 'k')
axis tight
hold on
x = xlim;
plot(x, [iis.thr(2) iis.thr(2)], '--r')
scatter(iis.peakPos(idx5) / fs / 60,...
    iis.peakPower(idx5), '*');
bsstamps = RestrictInts(bs.stamps, [idx2(1) - 20 * fs idx2(end) + 20 * fs]);
ylim(YY);
% yticks(Y);
if ~isempty(bsstamps)
    fill([bsstamps fliplr(bsstamps)] / fs / 60, [YY(1) YY(1) YY(2) YY(2)],...
        'k', 'FaceAlpha', 0.25,  'EdgeAlpha', 0);
end
% ylabel('Voltage [mV]')
xlabel('Time [m]')
xticks(idx3)
xticklabels([17.5 19.5])
xlim(idx3)
set(gca, 'TickLength', [0 0], 'YColor', 'none')
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bs identification

subplot(6, 2, [9, 11])
gscatter(varmat(:, 1), varmat(:, 3), gi, 'rk', '.', 1);
axis tight
hold on
gmPDF = @(x1, x2)reshape(pdf(gm, [x1(:) x2(:)]), size(x1));
ax = gca;
fcontour(gmPDF, [ax.XLim ax.YLim], 'HandleVisibility', 'off')
xlabel(['std log10(\sigma)]'])
ylabel('PC1')
legend({'Burst', 'Suppression'}, 'Location', 'northwest')
set(gca, 'TickLength', [0 0])
box off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aneStats
gidx = [ones(1, 30), ones(1, 30) * 2,...
    ones(1, 30) * 3, ones(1, 30) * 4];
c = [1 0 0; 1 0 1; 0 0 1; 0 1 1];
c2 = 'rmbc';

subplot(6, 2, [10, 12])
boxplot(as.deepFraction, gidx, 'PlotStyle', 'traditional',...
    'BoxStyle', 'outline', 'Color', c2, 'notch', 'off')
hold on
gscatter(gidx, [as.deepFraction(:)], gidx, c2)
legend off
xlabel('')
xticklabels({'WT', 'APPPS1', 'APPKi', 'FADx5'})
ylabel('Deep anesthesia duration [% total]')
box off
set(gca, 'TickLength', [0 0])
sigstar({[1, 3], [1, 4], [2, 4], [2, 3]}, [0.05, 0.01, 0.01, 0.05], 1);

figname = fullfile('F:\Daniel IIS 11.6.20\IIS\graphics', 'figS9');
export_fig(figname, '-tif', '-transparent', '-r300')
savePdf('figS9', 'F:\Daniel IIS 11.6.20\IIS', fh)

