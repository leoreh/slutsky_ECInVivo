function ripp = getRipples(varargin)

% Detects ripples from lfp: rectifies (squares) the signal, applies a
% moving average, standardizes, finds crossings of a threshold
%(in z-scores) and converts them to bouts by applying duration criterions.
% Alas, discards ripples with a low peak power. After detection, calculates
% various stats and plots a summary. based in part on
% bz_FindRipples but.
% 
% INPUT:
%   basepath        path to recording {pwd}
%   sig             numeric vec of lfp data. if empty will be loaded from
%                   basename.lfp according to rippCh.  
%   emg             numeric of of emg data. must be the same sampling
%                   frequency as sig. if empty will be loaded from
%                   basename.lfp according to emgCh.
%   fs          	numeric. sampling frequency of lfp file / data. 
%                   if empty will be extracted from session info (ce
%                   format)
%   rippCh          numeric vec. channels to load and average from lfp
%                   file {1}. if empty will be selected best on the ratio
%                   of mean to median within the passband. 
%   emgCh           numeric vec. emg channel in basename.lfp file. 
%   recWin          numeric 2 element vector. time to analyse in recording
%                   [s]. start of recording is marked by 0. {[0 Inf]}
%   graphics        logical. plot graphics {true} or not (false)
%   saveVar         logical. save variables (update spikes and save su)
%
% OUTPUT:
%   ripp            struct
%
% DEPENDENCIES:
%   binary2bouts
%   lfpFilter
%   Sync (buzcode)
%   SyncMap (buzcode)
%   PlotColorMap (buzcode)
%
% TO DO LIST:
%   finish graphics (done)
%   stats (done)
%   rate (done)
%   exclusion by emg noise
%   exclusion by spiking activity
%   exlude active periods when normalizing signal
%   allow user to input sig directly instead of loading from binary (done)
%   improve routine to select best ripple channel automatically
%
% 02 dec 21 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'sig', [], @isnumeric);
addOptional(p, 'emg', [], @isnumeric);
addOptional(p, 'recWin', [0 Inf]);
addOptional(p, 'rippCh', 1, @isnumeric);
addOptional(p, 'emgCh', [], @isnumeric);
addOptional(p, 'fs', [], @isnumeric);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
basepath    = p.Results.basepath;
sig         = p.Results.sig;
emg         = p.Results.emg;
recWin      = p.Results.recWin;
rippCh      = p.Results.rippCh;
emgCh       = p.Results.emgCh;
fs          = p.Results.fs;
graphics    = p.Results.graphics;
saveVar     = p.Results.saveVar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = pwd;
win = ones(11, 1) / 11;     % for moving average
shift = (length(win) - 1) / 2;
thr = [2 7];                % threshold of stds above the sig
limDur = [20, 150, 30];     % min, max, and inter dur limits for ripples [ms]
passband = [130 200];
binsizeRate = 30;           % binsize for calculating ripple rate [s]

% load session info
[~, basename] = fileparts(basepath);
sessionName = [basename, '.session.mat'];
if ~exist(sessionName, 'file')
    session = CE_sessionTemplate(pwd, 'viaGUI', false,...
        'force', true, 'saveVar', true);
else
    load(sessionName)
end

spkgrp = session.extracellular.spikeGroups.channels;
spkch = sort([spkgrp{:}]);
nchans = session.extracellular.nChannels;
if isempty(fs)
    fs = session.extracellular.srLfp;
end

fprintf('\ngetting ripples for %s\n', basename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(sig)    
    if isempty(rippCh)
        % if rippCh not specified, load 4 hr of data and find ch with best
        % SNR for ripples. note this subrutine occasionally points to the
        % channel with the highest movement artifacts / spikes. it is also
        % very time consuming. better to select manually.
        fprintf('selecting best ripple channel...\n')

        recDur = min([diff(recWin), 4 * 60 * 60]);
        sig = double(bz_LoadBinary([basename, '.lfp'], 'duration', recDur,...
            'frequency', fs, 'nchannels', nchans, 'start', recWin(1),...
            'channels', spkch, 'downsample', 1));
        
        sig_filt = filterLFP(sig, 'fs', fs, 'type', 'butter', 'dataOnly', true,...
            'order', 3, 'passband', passband, 'graphics', false);
        
        pow = fastrms(sig_filt, 15);
        ch_rippPowRatio = mean(pow) ./ median(pow);
        [~, rippCh] = max(ch_rippPowRatio);
    end
    
    % load all data and average across channels
    fprintf('loading lfp data...\n')
    sig = double(bz_LoadBinary([basename, '.lfp'], 'duration', diff(recWin),...
        'frequency', fs, 'nchannels', nchans, 'start', recWin(1),...
        'channels', rippCh, 'downsample', 1));
    if length(rippCh) > 1
        sig = mean(sig, 2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find and exclude bouts of high movement / active wake
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ignore_idx = false(length(sig), 1);
mov_idx = false(length(sig), 1);
wake_idx = false(length(sig), 1);

% -------------------------------------------------------------------------
% ALT 1: load emg data and find bouts of high activity
% load emg data from basename.lfp
if isempty(emg) && ~isempty(emgCh)
    emg = double(bz_LoadBinary([basename, '.lfp'], 'duration', diff(recWin),...
        'frequency', fs, 'nchannels', nchans, 'start', recWin(1),...
        'channels', emgCh, 'downsample', 1));
end

% make sure emg and sig (lfp) are the same length. if not, assumes
% discrepancy is due to the non-integer sampling frequency of tdt and
% fixes this by interpolating the shorter of the two signals.
lenDiff = length(emg) - length(sig);
if lenDiff ~= 0
    fprintf('emg and sig differ by %d samples (%.2f%%). interpolating...\n',...
        lenDiff, abs(lenDiff / length(sig) * 100))
    if lenDiff < 0
        tstamps_sig = [1 : length(sig)] / fs;
        emg = [interp1([1 : length(emg)] / fs, emg, tstamps_sig, 'pchip')]';        
    elseif lenDiff > 0
        tstamps_sig = [1 : length(emg)] / fs;
        sig = [interp1([1 : length(sig)] / fs, sig, tstamps_sig, 'pchip')]';
    end
end

% find bouts of high activity
if ~isempty(emg)
    fprintf('finding bouts of high activity...\n')   
    emg_rms = fastrms(emg, 15);
    mov_idx = emg_rms > prctile(emg_rms, 80);    
end

% ALT 2: use active wake from sleep states
ssfile = fullfile(basepath, [basename '.AccuSleep_states.mat']);
if exist(ssfile)    
    load(ssfile, 'ss')    
    [~, cfg_names, ~] = as_loadConfig([]);
    wake_stateIdx = find(strcmp(cfg_names, 'WAKE'));
    

    wake_inInt = InIntervals(ss.boutTimes{wake_stateIdx}, recWin);
    wake_bouts = ss.boutTimes{wake_stateIdx}(wake_inInt, :) * fs;  
    for iwake = 1 : size(wake_bouts, 1)
        wake_idx(wake_bouts(iwake, 1) : wake_bouts(iwake, 2)) = true;
    end
end
ignore_idx = wake_idx | mov_idx;

fprintf('preparing signal...\n')

% filter lfp data in ripple band
sig_filt = filterLFP(sig, 'fs', fs, 'type', 'butter', 'dataOnly', true,...
    'order', 3, 'passband', passband, 'graphics', false);

% remove samples of high movement 
nss = sig_filt;
nss(ignore_idx) = nan;

% sqaure, moving average and correct shift
nss = nss .^ 2;
% [nss, z] = filter(win, 1, nss);
% nss = [nss(shift + 1 : end, :); z(1 : shift, :)];

% standardize
nss = (nss - mean(nss, 'omitnan')) / std(nss, 'omitnan'); 

% debugging of mov_idx ----------------------------------------------------
debugFlag = 0;
if debugFlag
    fh = figure;
    sb1 = subplot(2, 1, 1);
    plot([1 : length(nss)] / 1250 / 60 / 60, nss);
    hold on
    scatter(find(nss > thr(1)) / 1250 / 60 / 60, 10 * ones(1, sum(nss > thr(1))), '.');
    xlabel('Time [h]')
    ylabel('NSS')
    legend({'', 'nss > thr'})
    
    sb2 = subplot(2, 1, 2);
    plot([1 : length(nss)] / 1250 / 60 / 60, emg_rms)
    hold on
    scatter(find(mov_idx) / 1250 / 60 / 60, 10 * ones(1, sum(mov_idx)), '.');
    xlabel('Time [h]')
    ylabel('EMG RMS')
    legend({'', 'emg > median'})
    linkaxes([sb1, sb2], 'x')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find ripples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find bouts with power > low threshold. correct for durations
limDur = limDur / 1000 * fs;
bouts = binary2bouts('vec', nss > thr(1), 'minDur', limDur(1),...
    'maxDur', limDur(2), 'interDur', limDur(3));

% discard ripples with a peak power < high threshold
peakPowNorm = zeros(size(bouts, 1), 1);
for ibout = 1 : size(bouts, 1)
    peakPowNorm(ibout) = max(nss(bouts(ibout, 1) : bouts(ibout, 2)));
end
dicard_idx = peakPowNorm < thr(2);
bouts(dicard_idx, :) = [];
peakPowNorm(dicard_idx) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find negative peak position for each ripple
peakPos = zeros(size(bouts, 1), 1);
peakPow = zeros(size(bouts, 1), 1);
for ibout = 1 : size(bouts, 1)
    [peakPow(ibout), peakPos(ibout)] =...
        min(sig_filt(bouts(ibout, 1) : bouts(ibout, 2)));
    peakPos(ibout) = peakPos(ibout) + bouts(ibout, 1) - 1;
end

% convert idx to seconds
nbouts = size(bouts, 1);
peakPos = peakPos / fs;
bouts = bouts / fs;

% instantaneous phase and amplitude
h = hilbert(sig_filt);
sig_phase = angle(h);
sig_amp = abs(h);
sig_unwrapped = unwrap(sig_phase);

% instantaneous frequency
tstamps = [1 / fs : 1 / fs : length(sig) / fs];
tstamps = tstamps(:);
dt = diff(tstamps);
t_diff = tstamps(1 : end - 1) + dt / 2;
d0 = diff(medfilt1(sig_unwrapped, 12)) ./ dt;
d1 = interp1(t_diff, d0, tstamps(2 : end - 1, 1));
sig_freq = [d0(1); d1; d0(end)] / (2 * pi);

% arrange maps
ripp.maps.durWin = [-75 75] / 1000;
nbins = floor(fs * diff(ripp.maps.durWin) / 2) * 2 + 1; % must be odd
centerBin = ceil(nbins / 2);
[r, i] = Sync([tstamps sig_filt], peakPos, 'durations', ripp.maps.durWin);
ripp.maps.ripp = SyncMap(r, i, 'durations', ripp.maps.durWin,...
    'nbins', nbins, 'smooth', 0);
[f, i] = Sync([tstamps, sig_freq], peakPos, 'durations', ripp.maps.durWin);
ripp.maps.freq = SyncMap(f, i, 'durations', ripp.maps.durWin, 'nbins', nbins, 'smooth', 0);
[a, i] = Sync([tstamps sig_phase], peakPos, 'durations', ripp.maps.durWin);
ripp.maps.phase = SyncMap(a, i,'durations', ripp.maps.durWin, 'nbins', nbins, 'smooth', 0);
[p, i] = Sync([tstamps sig_amp], peakPos, 'durations', ripp.maps.durWin);
ripp.maps.amp = SyncMap(p, i,'durations', ripp.maps.durWin, 'nbins', nbins, 'smooth', 0);

% more stats
ripp.maxFreq = max(ripp.maps.freq, [], 2);
ripp.peakFreq = ripp.maps.freq(:, centerBin);
ripp.peakAmp = ripp.maps.amp(:, centerBin);
ripp.dur = bouts(:, 2) - bouts(:, 1);

% acg and correlations
[ripp.acg.data, ripp.acg.t] = CCG(peakPos,...
    ones(length(peakPos), 1), 'binSize', 0.01);
ripp.corr.amp_freq = corrcoef(ripp.peakAmp, ripp.peakFreq);
ripp.corr.dur_freq = corrcoef(ripp.dur, ripp.peakFreq);
ripp.corr.dur_amp = corrcoef(ripp.dur, ripp.peakAmp);

% rate
[ripp.rate.rate, ripp.rate.binedges, ripp.rate.tstamps] =...
    times2rate(peakPos, 'binsize', binsizeRate, 'winCalc', [0, Inf],...
    'c2r', true);

% relation to sleep states
ssfile = fullfile(basepath, [basename '.AccuSleep_states.mat']);
if exist(ssfile)    
    load(ssfile, 'ss')    
    [cfg_colors, ~, ~] = as_loadConfig([]);

    ripp.states.stateNames = ss.labelNames;
    nstates = length(ss.boutTimes);
    sstates = [1 : 5];  % selected states
    
    for istate = sstates
        boutIdx = InIntervals(ss.boutTimes{istate}, recWin);
        if ~isempty(ss.boutTimes{istate})
            % rate in states
            [ripp.states.rate{istate}, ripp.states.binedges{istate},...
                ripp.states.tstamps{istate}] =...
                times2rate(peakPos, 'binsize', binsizeRate,...
                'winCalc', ss.boutTimes{istate}(boutIdx, :), 'c2r', true);
            
            % idx of rippels in state
            ripp.states.idx{istate} =...
                InIntervals(peakPos, ss.boutTimes{istate}(boutIdx, :));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    fh = figure;
    durPlot = [-50 50] / 1000;
    x = durPlot(1) : diff(durPlot) / nbouts : durPlot(2);
    histBins = 200;
    
    % examples on raw and filtered data
    sb1 = subplot(4, 3, [1, 2]);    
    idx_recMargin = 10 * fs;     % [s]
    rippSelect = 500 : 550;
    rippCenter = round((rippSelect(end) - rippSelect(1)) / 2) + rippSelect(1);
    idx_rec = round([peakPos(525) * fs - idx_recMargin :...
        peakPos(525) * fs + idx_recMargin]);
    plot(idx_rec / fs, sig(idx_rec), 'k')
    hold on
    plot(idx_rec / fs, sig_filt(idx_rec), 'm')
    yLimit = ylim;
    plot([peakPos(rippSelect), peakPos(rippSelect)], yLimit, 'b')
    plot([bouts(rippSelect, 1), bouts(rippSelect, 1)], yLimit, 'g')
    plot([bouts(rippSelect, 2), bouts(rippSelect, 2)], yLimit, 'r')
    xlim([peakPos(rippCenter) - 0.5, peakPos(rippCenter) + 0.5])
    legend('Raw LFP', 'Filtered')
    xlabel('Time [s]')
    
    % examples of ripples (filtered) superimposed
    sb3 = subplot(4, 3, 3);
    ripp_idx = randperm(nbouts, 1000);
    plot(((1 : nbins)' - ceil(nbins / 2)) / nbins * diff(durPlot),...
        ripp.maps.ripp(ripp_idx, :)', 'k');
    xlabel('Time [s]')
    
    % rate
    sb4 = subplot(4, 3, [4, 5]);
    plot(ripp.rate.tstamps / 60 / 60, ripp.rate.rate, 'k')
    xlabel('Time [h]')
    ylabel('Ripple Rate [Hz]')
    if exist(ssfile)
        hold on
        for istate = sstates
            ph = plot(ripp.states.tstamps{istate} / 60 / 60,...
                ripp.states.rate{istate}, '.', 'MarkerSize', 5);
            ph.Color = cfg_colors{istate};
        end
        xlabel('Time [h]')
        ylabel('Ripple Rate [Hz]')
        legend(["Total", ripp.states.stateNames{sstates}]);
    end
    ylim([0 3])
    
    % percent rippels in state
    sb6 = subplot(4, 3, 6);
    if exist(ssfile)
        pie(sum(cell2nanmat(ripp.states.idx), 1, 'omitnan'), ones(1, length(sstates)))
        hold on
        ph = findobj(sb6, 'Type', 'Patch');
        set(ph, {'FaceColor'}, flipud(cfg_colors(sstates)'))
    end
      
    % frequency map
    sb7 = subplot(4, 3, 7);
    PlotColorMap(ripp.maps.freq, 1, 'bar','on', 'cutoffs', [100 250], 'x', x);
    ylabel('Ripple No.')
    subtitle('Frequency');   
    
    % amplitude map
    sb8 = subplot(4, 3, 8);
    PlotColorMap(ripp.maps.amp, 1, 'bar','on', 'x', x);
    ylabel('Ripple No.')
    subtitle('Amplitude');
    
    % ACG
    sb9 = subplot(4, 3, 9);
    plotCCG(ripp.acg.data, ripp.acg.t);
    xlabel('Time [ms]')
    ylabel('Rate')
    
    % distribution of peak frequency
    sb10 = subplot(4, 3, 10);
    h = histogram(ripp.peakFreq, histBins, 'Normalization', 'probability');
    h.FaceColor = 'k';
    h.EdgeColor = 'none';
    xlabel('Peak Frequency [Hz]')
    ylabel('Probability')
    
    % distribution of peak amplitude
    sb11 = subplot(4, 3, 11);
    h = histogram(ripp.peakAmp, histBins, 'Normalization', 'probability');
    h.FaceColor = 'k';
    h.EdgeColor = 'none';
    xlabel('Peak Amp')
    ylabel('Probability')
    
    % distribution of ripple duration
    sb12 = subplot(4, 3, 12);
    h = histogram(ripp.dur * 1000, histBins, 'Normalization', 'probability');
    h.FaceColor = 'k';
    h.EdgeColor = 'none';
    xlabel('Ripple Duration [ms]')
    ylabel('Probability')
    
    % save figure
    figpath = fullfile(basepath, 'graphics');
    mkdir(figpath)
    figname = fullfile(figpath, sprintf('%s_ripples', basename));
    export_fig(figname, '-tif', '-transparent', '-r300')   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finalize and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ripp.info.rippCh = rippCh;
ripp.info.limDur = limDur;
ripp.info.recWin = recWin;
ripp.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
ripp.info.thr = thr;
ripp.bouts         = bouts + recWin(1);  
ripp.peakPos        = peakPos;
ripp.peakPow        = peakPow;
ripp.peakPowNorm    = peakPowNorm;

if saveVar      
    save([basepath, filesep, basename, '.ripp.mat'], 'ripp')
    
    % create ns files for visualization with neuroscope
    fs_dat = session.extracellular.sr;
    
    res = [ripp.bouts(:, 1); ripp.bouts(:, 2); ripp.peakPos] * fs_dat;
    [res, sort_idx] = sort(round(res));
    fid = fopen([basename, '.ripp.res.1'], 'w');
    fprintf(fid, '%d\n', res);
    rc = fclose(fid);
   
    clu = [ones(nbouts, 1); ones(nbouts, 1) * 2; ones(nbouts, 1) * 3];
    clu = clu(sort_idx);
    fid = fopen([basename, '.ripp.clu.1'], 'w');
    fprintf(fid, '%d\n', 3);
    fprintf(fid, '%d\n', clu);
    rc = fclose(fid);
end

end

% EOF