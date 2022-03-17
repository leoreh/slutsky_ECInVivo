
% spklfp_wrapper

% to do list
% adapt wavelet filtering (bz_WaveSpec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = pwd;
[~, basename] = fileparts(basepath);
cd(basepath)

% files
spksfile = fullfile(basepath, [basename, '.spikes.cellinfo.mat']);
unitsfile = fullfile(basepath, [basename, '.units.mat']);
lfpfile = fullfile(basepath, [basename, '.lfp']);
sessionfile = fullfile(basepath, [basename, '.session.mat']);
spklfpfile = fullfile(basepath, [basename, '.spklfp.mat']);
sleepfile = fullfile(basepath, [basename, '.sleep_states.mat']);

% load session info
if ~exist(sessionfile, 'file')
    session = CE_sessionTemplate(pwd, 'viaGUI', false,...
        'forceDef', false, 'forceL', false, 'saveVar', false);
else
    load(sessionfile)
end

% params from session info
nchans = session.extracellular.nChannels;
fs = session.extracellular.srLfp;

% load spikes
if exist(spksfile, 'file')
    load(spksfile, 'spikes')
end
nunits = length(spikes.times);

% load lfp
winCalc = [0 3 * 60 * 60];
recDur = diff(winCalc);
ch = 9 : 11;
sig = double(bz_LoadBinary(lfpfile, 'duration', diff(winCalc),...
    'frequency', fs, 'nchannels', nchans, 'start', winCalc(1),...
    'channels', ch, 'downsample', 1));
sig = mean(sig, 2);

% clip spike times
spktimes = cellfun(@(x) x(InIntervals(x, winCalc)),...
    spikes.times, 'uni', false);

% restrict to sleep states 
if exist(sleepfile, 'file')
    load(sleepfile, 'ss')
end
istate = 4;         % state of interest
spktimes = cellfun(@(x) x(InIntervals(x, ss.stateEpochs{istate})),...
    spktimes, 'uni', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single band analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frange = [0.5 2; 2, 4; 5, 11; 12, 18; 18, 30; 30, 50; 50, 80];

clear spklfp
for ifreq = 1 : size(frange, 1)
    
    fprintf('working on frequency %d of %d', ifreq, size(frange, 1))
    
    % filter lfp
    sig_filt = filterLFP(sig, 'fs', fs, 'type', 'butter', 'dataOnly', true,...
        'order', 3, 'passband', frange(ifreq, :), 'graphics', false);
    
    [spklfp(ifreq)] = spklfp_singleband('basepath', basepath, 'fs', fs,...
        'sig', sig_filt, 'spktimes', spktimes, 'frange', frange(ifreq, :),...
        'graphics', true, 'saveVar', true);    
end

% cat struct
sl = catfields(spklfp, 'catdef', 'addim');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spike triggered lfp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params
tstamps = [winCalc(1) : 1 / fs : winCalc(2) - 1 / fs];
durWin = [-500 500] / 1000;         % [sec]
nbinsMap = floor(fs * diff(durWin) / 2) * 2 + 1; % must be odd
centerBin = ceil(nbinsMap / 2);

% loop through cells
maxnspks = 1000;
for iunit = 1 : length(spktimes)
    
    bz_Counter(iunit, length(spktimes), 'spk-triggered LFP')

    % select spikes
    nspks = min(maxnspks, length(spktimes{iunit}));
    if nspks < 500
        sl.lfpmap{iunit} = nan;
        continue
    end
    spkidx = floor(linspace(1, length(spktimes{iunit}), maxnspks));

    % extract 1 sec of lfp sorrounding each cell
    [r, i] = Sync([tstamps' sig], spktimes{iunit}(spkidx), 'durations', durWin);
    sl.lfpmap{iunit} = SyncMap(r, i, 'durations', durWin,...
        'nbins', nbinsMap, 'smooth', 0);
end

% cat to single matrix
lfpmap = [];
for iunit = 1 : length(spktimes)
    if ~isnan(sl.lfpmap{iunit})
        lfpmap = [lfpmap; sl.lfpmap{iunit}];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save
save(fullfile(basepath, [basename, '.spklfp.mat']), 'sl')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wide band graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% center frequencies 
freq = sl.info.frange(1, :) + diff(sl.info.frange) / 2;

posnegcolor = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);

fh = figure;
th = tiledlayout(2, 3);

% syn mag correlation
nexttile
hold on
plot(log2(freq), sl.pop.synmag_r)
plot(log2(freq([1 end])), [0 0], 'k--')
box off
LogScale('x', 2)
axis tight
title('Syn - Mag Correlation')
xlabel('Frequency [Hz]');
ylabel('Rho')

% MRL per population
nexttile
plot(log2(freq), sl.pop.phase.mrl)
LogScale('x', 2)
axis tight
box off
xlabel('Frequency [Hz]');
ylabel('Mean Resultant Length')
legend({'RS', 'FS'})
title('Syn Phase')

% cell map of spkcount-mag correlation
nexttile
imagesc(1 : length(freq), 1 : nunits, sl.ratemag.r)
colormap(gca, posnegcolor)
ColorbarWithAxis([-0.2 0.2], 'rho')
xticks(1 : length(freq))
xticklabels(string(floor(freq)))
xlabel('Frequency [Hz]');
ylabel('Cell')
axis tight
axis xy
title('Spike - Mag Correlation')

% cell map of spike rate during in phase bins across frequencies
nexttile
data = squeeze(mean(mean(sl.ratemap.rate, 1, 'omitnan'), 3, 'omitnan'));
imagesc(sl.ratemap.phase_bins(:, 1), 1 : length(freq), data')
yticks(1 : length(freq))
yticklabels(string(floor(freq)))
ylabel('Frequency [Hz]')
xlabel('Phase [rad]')
xticks([0 : pi : 2 * pi])
xticklabels(string(0 : 2) + "pi")
title('Mean Rate Across Cells')

% mrl of all cells vs. frequency
nexttile
plot_boxMean(sl.phase.mrl, 'allPnts', true)
xticklabels(string(floor(freq)))
ylim([0 1])
xlabel('Frequency [Hz]')
ylabel('Mean Resultant Length')
title('MRL vs. Frequency')

% spike triggered lfp of 3 cells
xval = linspace(-0.5, 0.5, nbinsMap);
yval = sl.lfpmap{6};
for 
plot(xval, mean(yval, 'omitnan'));
hold on 
ylabel('LFP [mV]')
xlabel('Time [s]')
title('Spike Triggered LFP')


% save
figpath = fullfile(basepath, 'graphics', 'spklfp');
mkdir(figpath)
figname = fullfile(figpath, sprintf('%.spk-lfp_wideband', basename));
export_fig(figname, '-tif', '-transparent', '-r300')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bzcode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data
basepath = pwd;
[~, basename] = fileparts(basepath);
cd(basepath)

load([basename, '.spktimes.mat'])
load([basename, '.spikes.cellinfo.mat'])
load([basename, '.units.mat'])
load([basename, '.cell_metrics.cellinfo.mat'])

fs = 1250;

% clip spikes times to reasonable window
winCalc = [0, 3 * 60 * 60];
recDur = diff(winCalc);

% load lfp
sig = double(bz_LoadBinary([basename, '.lfp'], 'duration', recDur,...
    'frequency', fs, 'nchannels', 19, 'start', winCalc(1),...
    'channels', [9 : 11], 'downsample', 1));
sig = mean(sig, 2);

% filter
% sig_filt = filterLFP(sig, 'fs', fsLfp, 'type', 'butter', 'dataOnly', true,...
%     'order', 3, 'passband', [0.5 8], 'graphics', false);

% buzcode
subpops = cell(1, length(units.rs));
for iunit = 1 : length(subpops)
    if units.fs(iunit)
        subpops{iunit} = 'fs';
    elseif units.rs(iunit)
        subpops{iunit} = 'rs';
    end
end

lfp.data = sig;
lfp.timestamps = [winCalc(1) : 1 / fs : winCalc(2) - 1 / fs];
lfp.samplingRate = fs;
lfp.channels = 1;
[SpikeLFPCoupling] = bz_GenSpikeLFPCoupling(spikes, lfp,...
    'spikeLim', 500000, 'frange', [1 100], 'nfreqs', [20],...
    'cellclass', subpops, 'sorttype', 'rate', 'saveMat', true);


% -----------------------------
flfp.timestamps = [winCalc(1) : 1 / fs : winCalc(2) - 1 / fs];
sig_z = hilbert(sig);
flfp.amp = abs(sig_z);
flfp.phase = angle(sig_z);
flfp.samplingRate = fs;
bz_PowerPhaseRatemap(spikes, flfp)


% stimes = cellfun(@(x) x(InIntervals(x, winCalc))', spikes.times, 'uni', false);
% 
% % sort by fr
% [~, sidx] = sort(cellfun(@length, stimes, 'uni', true), 'descend');
% stimes = stimes(sidx);
% 
% % choose pyr
% stimes = stimes(units.idx(1, :));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wide band analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frange = [1 100];
nfreq = 20;
[spklfp] = spklfp_wideband('basepath', basepath, 'fs', fs,...
    'sig', sig, 'spktimes', spktimes, 'frange', frange,...
    'nfreq', nfreq, 'srtUnits', true, 'graphics', true,...
    'saveVar', true);

