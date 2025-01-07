function s = spklfp_wrapper(varargin)

% INPUT
%   basepath    char. fullpath to recording folder {pwd}
%   winCalc     cell array of n x 2 mats of intervals.
%               metrices will be calculated for each cell by limiting
%               spktimes to the intervals. can be for example
%               ss.boutTimes. must be the same units as spikes.times
%               (e.g. [s])
%   ch          numeric. lfp channels to load. can be a vector for averaging
%   frange      2 x n mat of frequency bands to analyze
%   saveVar     logical {true}
%   graphics    logical {true}
%
% CALLS
%   spklfp_singleband
%
% TO DO LIST
%   # adapt wavelet filtering (bz_WaveSpec)
%
% 25 feb 22 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar)
addParameter(p, 'winCalc', {[0 Inf]}, @iscell)
addParameter(p, 'ch', 1, @isnumeric)
addParameter(p, 'frange', [], @isnumeric)
addParameter(p, 'graphics', true, @islogical)
addParameter(p, 'saveVar', true, @islogical)

parse(p, varargin{:})
basepath        = p.Results.basepath;
winCalc         = p.Results.winCalc;
ch              = p.Results.ch;
frange          = p.Results.frange;
graphics        = p.Results.graphics;
saveVar         = p.Results.saveVar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% files and params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% frequency bands of interest
if isempty(frange)
    frange = [0.5 2; 2, 4; 5, 11; 12, 18; 18, 30; 30, 50; 50, 80];
end

% center frequencies
freq = frange(:, 1) + [diff(frange') / 2]';

% files
[~, basename] = fileparts(basepath);
cd(basepath)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iwin = 1 : length(winCalc)
    
    % load lfp
    winstart = min(winCalc{iwin});
    winend = max(winCalc{iwin});
    recDur = winend - winstart;
    sig = double(bz_LoadBinary(lfpfile, 'duration', recDur,...
        'frequency', fs, 'nchannels', nchans, 'start', winstart,...
        'channels', ch, 'downsample', 1));
    sig = mean(sig, 2);

    % clip spike times
    spktimes = cellfun(@(x) x(InIntervals(x, winCalc{iwin})),...
        spikes.times, 'uni', false);
    
    % adjust spktimes so that they correspond to signal length
    spktimes = cellfun(@(x) x - winCalc{iwin}(1),...
        spktimes, 'uni', false);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % single band analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clear spklfp
    for ifreq = 1 : size(frange, 1)

        fprintf('working on frequency %d of %d', ifreq, size(frange, 1))

        % filter lfp
        sig_filt = filterLFP(sig, 'fs', fs, 'type', 'butter', 'dataOnly', true,...
            'order', 3, 'passband', frange(ifreq, :), 'graphics', false);

        [spklfp(ifreq)] = spklfp_singleband('basepath', basepath, 'fs', fs,...
            'sig', sig_filt, 'spktimes', spktimes, 'frange', frange(ifreq, :),...
            'graphics', false, 'saveVar', false);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % spike triggered lfp
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    stlname = sprintf('stl_wideband_win%d', iwin);
    stl = spklfp_stl('basepath', basepath, 'fs', fs,...
            'sig', sig, 'spktimes', spktimes, 'mapWin', [],...
            'graphics', true, 'saveVar', stlname);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cat across frequencies
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    s(iwin) = catfields(spklfp, 'catdef', 'addim');
    s(iwin).info.winCalc = winCalc{iwin};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % wide band graphics
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    setMatlabGraphics(false)
    fh = figure;
    th = tiledlayout(2, 3, 'TileSpacing', 'Compact');
    figname = sprintf('%s; %d-%d hr',...
        basename, round(winCalc{iwin} / 60 / 60));
    title(th, figname, 'Interpreter', 'none')
    posnegcolor = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);

    % syn mag correlation
    nexttile
    hold on
    plot(log2(freq), s(iwin).pop.synmag_r)
    plot(log2(freq([1 end])), [0 0], 'k--')
    box off
    LogScale('x', 2)
    axis tight
    title('Syn - Mag Correlation')
    xlabel('Frequency [Hz]');
    ylabel('Rho')

    % MRL per population
    nexttile
    plot(log2(freq), s(iwin).pop.phase.mrl)
    LogScale('x', 2)
    axis tight
    box off
    xlabel('Frequency [Hz]');
    ylabel('Mean Resultant Length')
    legend({'RS', 'FS'})
    title('Syn Phase')

    % cell map of spkcount-mag correlation
    nexttile
    imagesc(1 : length(freq), 1 : nunits, s(iwin).ratemag.r)
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
    data = squeeze(mean(mean(s(iwin).ratemap.rate, 1, 'omitnan'), 3, 'omitnan'));
    imagesc(s(iwin).ratemap.phase_bins(:, 1), 1 : length(freq), data')
    yticks(1 : length(freq))
    yticklabels(string(floor(freq)))
    ylabel('Frequency [Hz]')
    xlabel('Phase [rad]')
    xticks([0 : pi : 2 * pi])
    xticklabels(string(0 : 2) + "pi")
    title('Mean Rate Across Cells')

    % mrl of all cells vs. frequency
    nexttile
    plot_boxMean(s(iwin).phase.mrl, 'allPnts', true)
    xticklabels(string(floor(freq)))
    ylim([0 1])
    xlabel('Frequency [Hz]')
    ylabel('Mean Resultant Length')
    title('MRL vs. Frequency')

    % save figure
    figpath = fullfile(basepath, 'graphics', 'spklfp');
    mkdir(figpath)
    figname = sprintf('%s_spklfp_wideband_win%d', basename, iwin);
    figname = fullfile(figpath, figname);
    export_fig(figname, '-tif', '-transparent', '-r300')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save var
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save
save(fullfile(basepath, [basename, '.spklfp.mat']), 's', '-v7.3')


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bzcode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % data
% basepath = pwd;
% [~, basename] = fileparts(basepath);
% cd(basepath)
%
% load([basename, '.spktimes.mat'])
% load([basename, '.spikes.cellinfo.mat'])
% load([basename, '.units.mat'])
% load([basename, '.cell_metrics.cellinfo.mat'])
%
% fs = 1250;
%
% % clip spikes times to reasonable window
% winCalc = [0, 3 * 60 * 60];
% recDur = diff(winCalc);
%
% % load lfp
% sig = double(bz_LoadBinary([basename, '.lfp'], 'duration', recDur,...
%     'frequency', fs, 'nchannels', 19, 'start', winCalc(1),...
%     'channels', [9 : 11], 'downsample', 1));
% sig = mean(sig, 2);
%
% % filter
% % sig_filt = filterLFP(sig, 'fs', fsLfp, 'type', 'butter', 'dataOnly', true,...
% %     'order', 3, 'passband', [0.5 8], 'graphics', false);
%
% % buzcode
% subpops = cell(1, length(units.rs));
% for iunit = 1 : length(subpops)
%     if units.fs(iunit)
%         subpops{iunit} = 'fs';
%     elseif units.rs(iunit)
%         subpops{iunit} = 'rs';
%     end
% end
%
% lfp.data = sig;
% lfp.timestamps = [winCalc(1) : 1 / fs : winCalc(2) - 1 / fs];
% lfp.samplingRate = fs;
% lfp.channels = 1;
% [SpikeLFPCoupling] = bz_GenSpikeLFPCoupling(spikes, lfp,...
%     'spikeLim', 500000, 'frange', [1 100], 'nfreqs', [20],...
%     'cellclass', subpops, 'sorttype', 'rate', 'saveMat', true);
%
%
% % -----------------------------
% flfp.timestamps = [winCalc(1) : 1 / fs : winCalc(2) - 1 / fs];
% sig_z = hilbert(sig);
% flfp.amp = abs(sig_z);
% flfp.phase = angle(sig_z);
% flfp.samplingRate = fs;
% bz_PowerPhaseRatemap(spikes, flfp)
%
%
% % stimes = cellfun(@(x) x(InIntervals(x, winCalc))', spikes.times, 'uni', false);
% %
% % % sort by fr
% % [~, sidx] = sort(cellfun(@length, stimes, 'uni', true), 'descend');
% % stimes = stimes(sidx);
% %
% % % choose pyr
% % stimes = stimes(units.idx(1, :));
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % wide band analysis
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% frange = [1 100];
% nfreq = 20;
% [spklfp] = spklfp_wideband('basepath', basepath, 'fs', fs,...
%     'sig', sig, 'spktimes', spktimes, 'frange', frange,...
%     'nfreq', nfreq, 'srtUnits', true, 'graphics', true,...
%     'saveVar', true);
%
