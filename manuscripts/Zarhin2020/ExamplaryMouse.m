% select mouse
i = 1;      % grp
j = 4;     % mouse

% params
graphics = true;
loaddata = false;
savefig = false; 

basepath{1} = 'E:\Data\Others\DZ\IIS\WT';
basepath{2} = 'E:\Data\Others\DZ\IIS\APPPS1';
basepath{3} = 'E:\Data\Others\DZ\IIS\APPKi';
basepath{4} = 'E:\Data\Others\DZ\IIS\FADx5';
cd(basepath{i})
filename = dir('*.abf');
files = natsort({filename.name});
nfiles = 1 : length(files);
[~, grpname] = fileparts(basepath{i});
[~, basename] = fileparts(files{nfiles(j)});

fs = 1250;
binsize = (2 ^ nextpow2(5 * fs));
smf = 7;
marg = 0.05;
thr = [0 0];
ch = 1;


if loaddata
% load lfp
lfp = getLFP('basepath', basepath{i}, 'ch', ch, 'chavg', {},...
    'fs', 1250, 'interval', [0 inf], 'extension', 'abf', 'pli', true,...
    'savevar', false, 'force', false, 'basename', basename);
sig = double(lfp.data(:, ch));

% load iis
iis = getIIS('sig', sig, 'fs', fs, 'basepath', basepath{i},...
    'graphics', false, 'saveVar', false, 'binsize', binsize,...
    'marg', marg, 'basename', basename, 'thr', thr, 'smf', 7,...
    'saveFig', false, 'forceA', false, 'spkw', false, 'vis', false);
wvstamps = linspace(-marg, marg, floor(marg * fs) * 2 + 1);

% load bs
vars = {'std', 'max', 'sum'};
bs = getBS('sig', sig, 'fs', fs, 'basepath', basepath{i}, 'graphics', false,...
    'saveVar', false, 'binsize', 0.5, 'BSRbinsize', binsize, 'smf', smf,...
    'clustmet', 'gmm', 'vars', vars, 'basename', basename,...
    'saveFig', false, 'forceA', false, 'vis', false);

% load delta
load([grpname '_as.mat'])
load([basename '.ep.mat'])
end

% idx for zoomin 
minmarg = 1.5;
midsig = 43;
idx = round((midsig - minmarg) * fs * 60 : (midsig + minmarg) * fs * 60);

if graphics    
    fh = figure('Visible', 'on');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    suptitle(basename)
    
    % spectrogram
    sb1 = subplot(4, 2, 1 : 2);
    specBand('sig', sig, 'graphics', true, 'binsize', binsize, 'smf', smf, 'normband', true);
    set(gca, 'TickLength', [0 0], 'XTickLabel', [],...
        'Color', 'none', 'XColor', 'none')
    ylim([0 100])
    title('')
    
    % bsr
    sb2 = subplot(4, 2, 3 : 4);
    plot(bs.cents / fs / 60, bs.bsr, 'k', 'LineWidth', 1)
    hold on
    plot(bs.cents / fs / 60, ep.dband, 'b', 'LineWidth', 1)
    ylim([0 1])
    Y = ylim;
    set(gca, 'TickLength', [0 0], 'YTickLabel', [], 'XTickLabel', [],...
        'Color', 'none', 'XColor', 'none', 'YColor', 'none')
    box off
    legend({'BSR', 'Delta'})
    if ~isempty(ep.deep_stamps)
        fill([ep.deep_stamps fliplr(ep.deep_stamps)]' / fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
            'b', 'FaceAlpha', 0.2,  'EdgeAlpha', 0, 'HandleVisibility', 'off');
    end
    if ~isempty(ep.sur_stamps)
        fill([ep.sur_stamps fliplr(ep.sur_stamps)]' / fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
            'g', 'FaceAlpha', 0.2,  'EdgeAlpha', 0, 'HandleVisibility', 'off');
    end
    axis tight

    % raw
    sb3 = subplot(4, 2, 5 : 6);
    plot(lfp.timestamps / 60, sig, 'k', 'LineWidth', 1)
    hold on
    set(gca, 'TickLength', [0 0], 'YTickLabel', [],...
        'Color', 'none')
    box off
    ylabel('LFP [mV]')
    Y = ylim;
    idx3 = [idx(1) idx(end)] / fs / 60;
    fill([idx3 fliplr(idx3)]', [Y(1) Y(1) Y(2) Y(2)],...
        'r', 'FaceAlpha', 0.2,  'EdgeAlpha', 0, 'HandleVisibility', 'off');
        axis tight
    xlabel('Time [m]')
    
    % zoom in
    subplot(4, 2, 8);
    idx2 = iis.peakPos > idx(1) & iis.peakPos < idx(end);
    plot(lfp.timestamps(idx) / 60, sig(idx), 'k')
    axis tight
    hold on
    x = xlim;
    plot(x, [iis.thr(2) iis.thr(2)], '--r')
    scatter(iis.peakPos(idx2) / fs / 60,...
        iis.peakPower(idx2), '*');
    bsstamps = intervals(bs.stamps).intersect(intervals([idx(1) idx(end)])).ints;
    Y = ylim;
    if ~isempty(bsstamps)
        fill([bsstamps fliplr(bsstamps)] / fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
            'k', 'FaceAlpha', 0.25,  'EdgeAlpha', 0);
    end
    ylabel('Voltage [mV]')
    xlabel('Time [m]')
    xticks(idx3)
    set(gca, 'TickLength', [0 0])
    box off
    title('IIS')    
    
    % iis waveforms
    subplot(4, 2, 7)
    plot(wvstamps * 1000, iis.wv)
    ylabel('Voltage [mV]')
    xlabel('Time [ms]')
    axis tight
    xticks([-marg, 0, marg] * 1000);
    set(gca, 'TickLength', [0 0])
    box off
    title('IIS waveform')
    
    % mean + std waveform
    axes('Position',[.1315 .11 .13 .08])
    box on
    stdshade(iis.wv, 0.5, 'k', wvstamps)
    axis tight
    set(gca, 'TickLength', [0 0], 'YTickLabel', [], 'XTickLabel', [],...
        'XColor', 'none', 'YColor', 'none', 'Color', 'none')
    title(sprintf('n = %d', size(iis.wv, 1)));
    box off
    
    linkaxes([sb1, sb2, sb3], 'x');

    if saveFig
        figname = [basename];
        export_fig(figname, '-tif', '-transparent')
        % savePdf(figname, basepath, ff)
    end
    
end   
    
    