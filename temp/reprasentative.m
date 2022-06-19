
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file
basepath = 'F:\Data\lh100\lh100_220403_100052';
[~, basename] = fileparts(basepath);

% load data
varsFile = ["fr"; "spikes"; "datInfo"; "session";...
    "units"; "sleep_states"];
varsName = ["fr"; "spikes"; "datInfo"; "session";...
    "units"; "ss"];
v = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
    'varsName', varsName);

% session params
nchans = v.session.extracellular.nChannels;
fs = v.session.extracellular.sr;

% sleep signals
filename = fullfile(basepath, [basename, '.sleep_sig.mat']);
load(filename, 'emg');
load(filename, 'eeg');
yLimit_emg = [prctile(emg, 0.05), prctile(emg, 99.95)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% re-calc stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% state params
cfg = as_loadConfig();
sstates = [6, 3, 5, 4, 2, 1];

% re-calc state epochs for better visualization
minDur = [10, 5, 5, 10, 5, 5];
minDur = 4;
interDur = 10;
labels = v.ss.labels;
labels(labels == 6) = 5;
labels(labels == 3) = 4;
[stateEpochs, ~] = as_epochs('labels', labels,...
    'minDur', minDur, 'interDur', interDur);

% subsample emg
fs_eeg = 250;
emg = emg(1 : 1250 / fs_eeg : end);
emg_tstamps = [1 : length(emg)] / fs_eeg;

% subsample eeg and re-calc spec
eeg = eeg(1 : 1250 / fs_eeg : end);
spec = calc_spec('sig', eeg, 'fs', fs_eeg, 'graphics', false, 'saveVar', false,...
    'padfft', -1, 'winstep', 1, 'logfreq', true, 'ftarget', [],...
    'ch', [{1}], 'force', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% windows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate 4 windows of different durations; 1 for plotting and 3 for
% x-limit (e.g., raster that is zoom in view of spectrogram). plotStart
% determines how many representatives to create (e.g. before and after
% injection)

winDur(1) = 60 * 60;        % for plotting and scrolling through the data
winDur(2) = 10 * 60;        % for hypnogram, spec, and emg
winDur(3) = 20;             % for raster plot
winDur(4) = 0.5;            % for raw data

% one row for every window. one column for different sections of the
% recording
winStart(1, 1) = 60 * 60;
winStart(2, 1) = [1150];        % relative to window 1
winStart(3, 1) = [1310];        % relative to window 1
winStart(4, 1) = [3.8];           % relative to window 3

winStart(1, 2) = v.session.general.timepnt + 130 * 60;
winStart(2, 2) = [1050];
winStart(3, 2) = [1140];
winStart(4, 2) = [1];           

nwin = size(winStart, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fig params
nrows = 5;
fh = figure;
th = tiledlayout(nrows, nwin);
th.TileSpacing = 'compact';
th.Padding = 'none';

% plot
clear ypush
for iwin = 1 : nwin
    
    % win 
    winLim = [winStart(:, iwin)'; winStart(:, iwin)' + winDur]';

    % hypnogram -----------------------------------------------------------
    axh(iwin) = nexttile(iwin);
    hold on
    for istate = 1 : length(sstates)
        winIdx = InIntervals(stateEpochs{sstates(istate)}, winLim(1, :));
        sepochs = stateEpochs{sstates(istate)}(winIdx, :);
        if ~isempty(sepochs)
            plot(sepochs' - winStart(1, iwin), 1 * ones(size(sepochs))',...
                'color', cfg.colors{sstates(istate)}, 'LineWidth', 25)
        end
    end
    set(gca, 'ytick', [])
    set(gca, 'YColor', 'none')
    pbaspect([10, 1, 1])

    % spectrogram
    axh(iwin + nwin * 1) = nexttile(iwin + nwin * 1);
    winIdx = InIntervals(spec.tstamps, winLim(1, :));
    specTmp = spec;
    specTmp.s = spec.s(winIdx, :);
    specTmp.tstamps = spec.tstamps(winIdx) - winStart(1, iwin);
    plot_spec(specTmp, 'ch', 1, 'logfreq', true, 'saveFig', false,...
        'axh', axh(iwin + nwin * 1), 'xtime', 1)
    ylim([0 60])
    yticks([1, 10, 50])

    % emg -----------------------------------------------------------------
    axh(iwin + nwin * 2) = nexttile(iwin + nwin * 2);
    hold on
    winIdx = InIntervals(emg_tstamps, winLim(1, :));
    plot(emg_tstamps(winIdx) - winStart(1, iwin), emg(winIdx))
    ylim(yLimit_emg)

    % add marking of narrow window
    yLimit = ylim;
    plot(winLim(3, :), [yLimit(2) / 2 yLimit(2) / 2], 'k', 'LineWidth', 2)

    % link hypnogram, spec, and emg axis
    linkaxes(axh(iwin : nwin : (nrows - 2) * nwin), 'x')
    xlim(winLim(2, :))

    % raster plot ---------------------------------------------------------
    % axis 
    axh(iwin + nwin * 3) = nexttile(iwin + nwin * 3);
    hold on
    
    % add marking of narrow window
    yLimit = ylim;
    plot(winLim(4, :) + winStart(3, iwin), [0, 0], 'k', 'LineWidth', 2)

    % arrange units
    [~, unitsIdx] = sortrows(v.units.clean', 'descend');
    nunits = sum(v.units.clean, 2);
    
    % plot all units in blue
    spktimes = cellfun(@(x) [x(InIntervals(x, winLim(1, :)))  - winStart(1, iwin)]',...
        v.spikes.times, 'uni', false)';
    spktimes = spktimes(unitsIdx(1 : sum(nunits)));   
    LineFormat.Color = [0.1 0.1 0.8];
    lineFormat.LineWidth = 2;
    plotSpikeRaster(spktimes,...
        'PlotType', 'vertline', 'LineFormat', LineFormat);
    
    % plot fs units in red    
    LineFormat.Color = [0.8 0.1 0.1];
    for iunit = 1 : nunits(1)
        spktimes{iunit} = 0;
    end  
    plotSpikeRaster(spktimes,...
        'PlotType', 'vertline', 'LineFormat', LineFormat);    
    
    % axis params
    ylabel('Unit #')
    xlim(winLim(3, :))

    % raw data ___---------------------------------------------------------
    % load raw data
    fs_raw = fs;
    raw = double(bz_LoadBinary([basename, '.dat'], 'duration', winDur(3),...
        'frequency', fs, 'nchannels', nchans, 'start', winStart(3, iwin),...
        'channels', [13 : 16], 'downsample', fs / fs_raw));
    raw_tstamps = [1 : size(raw, 1)] / fs_raw + winStart(3, iwin);
    
    % plot
    axh(iwin + nwin * 4) = nexttile(iwin + nwin * 4);
    if ~exist('ypush', 'var')
        ypush = median(prctile(raw, 0.1)) * [0 : size(raw, 2) - 1];
    end
    plot(raw_tstamps, raw + ypush, 'k')
    xlim(winLim(4, :) + winStart(3, iwin))

end

% save
set(gcf,'renderer','Painters')
figname = fullfile(basepath, 'graphics', [basename, '_rep']);
export_fig(figname, '-pdf', 'renderer', 'painters')