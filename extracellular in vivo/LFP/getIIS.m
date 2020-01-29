function [iis] = getIIS(varargin)

% detects inter-ictal spikes from LFP.
%
% INPUT
%   sig         signal for detection
%   fs          sampling frequency
%   smf         smooth factor for rate [bins] {6}.
%   binsize     scalar {60} in [s]. for rate calculation
%   marg        scalar {0.1} in [s]. time margin for clipping spikes
%   thr         vector of two elements. first is thr in [z-scores] and
%               second is thr in [mV]. if one specified than the other is
%               merely calculated. if both specificed than both used for
%               detection.
%   spkw        logical {false}. calculate max frequency for each spike.
%   basepath    recording session path {pwd}
%   basename    string. if empty extracted from basepath
%   graphics    logical {true}. plot figure
%   saveVar     logical {true}. save variable
%   saveFig     logical {true}. save figure
%   vis         logical. figure visible {1} or not
%   forceA      logical {false}. force analysis even if .mat exists
%   filtspk     logical {true}. filter spikes or not.
%
% OUTPUT
%   iis         struct with fields:
%       stamps      n x 2 mat where n is the number of events and the
%                   columns represent the start and end of each event [samps]
%       peaks       the timestamp for each peak [samps]
%       peakPower   the voltage at each peak [uV]
%       thr         10 z-scores in [mV]
%       edges       of bins for bsr
%       cents       of bins for bsr
%       <params>    as in input
%
% TO DO LIST
%       # investigate thr (done)
%       # filter IIS by removing average waveform (done)
%       # find width via wavelet and after hyperpolirazation (done)
%       # manage bidirectional thresholding
%       # remove spikes with more than one peak in window to produce a
%       clean mean waveform
%
% CALLS
%       calcFR
%
% 02 jan 20 LH      UPDATES:
% 13 jan 20 LH      filtspk
% 19 jan 20 LH      thr in mV and z-scores
% 22 jan 20 LH      max frequency at peakPos rather than throught spike
% 27 jan 20 LH      trough between peaks instead of refractory period
% 28 jan 20 LH      changed iis.rate to counts / bin instead of spikes / m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'sig', [], @isnumeric)
addParameter(p, 'fs', 1250, @isnumeric)
addParameter(p, 'smf', 6, @isnumeric)
addParameter(p, 'binsize', 300, @isnumeric)
addParameter(p, 'marg', 0.05, @isnumeric)
addParameter(p, 'thr', [10 0], @isnumeric)
addParameter(p, 'spkw', false, @islogical)
addParameter(p, 'basepath', pwd, @isstr);
addParameter(p, 'basename', [], @isstr);
addParameter(p, 'graphics', true, @islogical)
addParameter(p, 'saveVar', true, @islogical);
addParameter(p, 'saveFig', true, @islogical);
addParameter(p, 'vis', true, @islogical);
addParameter(p, 'forceA', false, @islogical);
addParameter(p, 'filtspk', true, @islogical);

parse(p, varargin{:})
sig = p.Results.sig;
fs = p.Results.fs;
smf = p.Results.smf;
marg = p.Results.marg;
thr = p.Results.thr;
spkw = p.Results.spkw;
binsize = p.Results.binsize;
basepath = p.Results.basepath;
basename = p.Results.basename;
graphics = p.Results.graphics;
saveVar = p.Results.saveVar;
saveFig = p.Results.saveFig;
vis = p.Results.vis;
forceA = p.Results.forceA;
filtspk = p.Results.filtspk;

% params
margs = floor(marg * fs);           % margs [samples]; marg [ms]
wvstamps = linspace(-marg, marg, margs * 2 + 1);
tstamps = [1 : length(sig)]' / fs;
interDur = round(fs * 0.025);       % samples
lowthr = 0.2;                       % mV

% initialize output
iis.edges = []; iis.cents = []; iis.peakPos = []; iis.peakPower = [];
iis.wv = []; iis.rate = []; iis.filtered = []; iis.fs = fs;
iis.binsize = binsize; iis.spkw = []; iis.maxf = []; iis.out = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if mat already exists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(basename)
    [~, basename] = fileparts(basepath);
end
filename = [basepath, '\', basename, '.iis.mat'];
if ~forceA
    if exist(filename)
        load(filename)
        fprintf('\n loading %s \n', filename)
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detect spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% threshold
if thr(1) == 0
    thresholded = sig > thr(2);
    % find thr in z scores
    thr(1) = (thr(2) - mean(sig)) / std(sig);
elseif thr(2) == 0
    thresholded = zscore(sig) > thr(1);
    % find thr in mV
    thr(2) = mean(sig) + thr(1) * std(sig);
else
    thresholded = zscore(sig) > thr(1) & sig > thr(2);
end
iis.thr = thr;
iie = find([0; diff(thresholded) > 0]);

% return if no IIS are detected
nspks = length(iie);
if isempty(iie)
    fprintf('\nno inter-ictal spikes found\n');
    iis.rate = zeros(1, floor(length(sig) / binsize));
    iis.cents = zeros(1, floor(length(sig) / binsize));
    if saveVar
        save(filename, 'iis')
    end
    return
end

% remove crossings that are too close together. this effectively limits the
% maximum IIS burst frequency to 1 / interDur. the division of interDur is
% so that the search for max later on extends beyond the deleted
% crossings
ii = find(diff(iie) < interDur / 2);
while ~isempty(ii)
    iie(ii + 1) = [];
    ii = find(diff(iie) < interDur / 2);
end

% select local maximum and clip spikes 
peak = zeros(length(iie), 1);
pos = zeros(length(iie), 1);
seg = zeros(length(iie), length(wvstamps));
rmidx = [];
for i = 1 : length(iie)
    % remove last / first iis if incomplete
    if iie(i) + margs > length(sig)
        rmidx = [rmidx, i];
        continue
    elseif iie(i) - margs < 1
        rmidx = [rmidx, i];
        continue
    end
    % adjust crossing to local maximum and then clip
    localseg = sig(iie(i) - interDur : iie(i) + interDur);
    [peak(i), pos(i)] = max(localseg);
    pos(i) = iie(i) - interDur + pos(i);
    seg(i, :) = sig(pos(i) - margs : pos(i) + margs);
end
% if there is no trough between peaks select highest peak. this is
% instead of demanding for a refractory period via interDur
for i = 1 : length(pos) - 1
    low = min(sig(pos(i) : pos(i + 1)));
    if low > lowthr
        rmidx = [rmidx; i + (peak(i) > peak(i + 1))];
    end
end
seg(rmidx, :) = [];
peak(rmidx) = [];
pos(rmidx) = [];

mwv = mean(seg, 1)';
nspks = length(pos);
fprintf('\n%d inter-ictal spikes \n', nspks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% power spectrum via wavelet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% average wavelet coefficients for each spike. note this produces very
% similar results to the coefficients obtained from the average spike
% waveform. individual coefficients may still be carried out if the width
% of each spike is to be calculated

% filter bank
fb = cwtfilterbank('SignalLength', size(seg, 2), 'VoicesPerOctave', 32,...
    'SamplingFrequency', fs, 'FrequencyLimits', [1 625]);
if spkw
    for i = 1 : size(seg, 1)
        % ALT 1: spike width by inverse of max wavelet
        [cfs(i, :, :), f, coi] = cwt(seg(i, :), 'FilterBank', fb);
        [~, ifreq] = max(abs(squeeze(cfs(i, :, :))), [], 1);
        iis.maxf(i) = f(ifreq(round(size(seg, 2) / 2)));
        iis.spkw(i) = 1000 / iis.maxf(i);
        
        % ALT 2: spike width at half maximum
        % [~, ~, w(i)] = findpeaks(-seg(i, :), fs, 'NPeaks', 1,...
        %    'WidthReference', 'halfheight', 'MinPeakHeight', abs(thr(2)));
    end
    cfs = squeeze(mean(abs(cfs), 1));
    
    % idx to suspecious spikes
    % iis.out = find(isoutlier(iis.spkw, 'ThresholdFactor', 2));
    iis.out = find(iis.spkw < 5);
else
    [cfs, f, coi] = cwt(mwv, 'FilterBank', fb);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[iis.rate, iis.edges, iis.cents] = calcFR(pos, 'winCalc', [1, length(sig)],...
    'binsize', binsize, 'smet', 'none', 'c2r', false);
% this is super dangerous because if binsize < 1 min than the rate will
% effectively be greater than the number of counts
% iis.rate = iis.rate * fs * 60;      % convert counts in bins to 1 / min
iis.rate = movmean(iis.rate, smf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ccg, tccg] = CCG({pos / fs}, [], 'duration', 10,...
    'binSize', 0.1);
[ccg2, tccg2] = CCG({pos / fs}, [], 'duration', 30,...
    'binSize', 0.3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filtered = sig;
if filtspk
    for i = 1 : size(seg, 1)
        idx = [pos(i) - margs : pos(i) + margs];
        filtered(idx) = filtered(idx) - mwv;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange struct output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iis.wv = seg;
iis.filtered = filtered;
iis.peakPos = pos;
iis.peakPower = peak;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveVar
    save(filename, 'iis')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    if vis
        fh = figure;
    else
        fh = figure('Visible', 'off');
    end
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    suptitle(basename)
    
    % raw and iis
    subplot(3, 4, 1 : 2)
    plot(tstamps / 60, sig ,'k')
    axis tight
    hold on
    plot(xlim, [thr(2) thr(2)], '--r')    
    ylabel('Voltage [mV]')
    yyaxis right
    plot(iis.cents / fs / 60, iis.rate, 'b', 'LineWidth', 3)
    ylabel('Rate [spikes / bin]')
    set(gca, 'TickLength', [0 0])
    box off
    title('Raw signal and IIS rate')
    
    % detection
    subplot(3, 4, 5 : 6)
    plot(tstamps / 60, zscore(sig), 'k')
    hold on
    axis tight
    plot(xlim, [thr(1) thr(1)], '--r')
    plot([pos pos] / fs / 60, [-10 -1], '--g', 'LineWidth', 2)
    xlabel('Time [m]')
    ylabel('Z-score')
    set(gca, 'TickLength', [0 0])
    box off
    title('IIS detection')
    
    % zoom in
    subplot(3, 4, 9 : 10)
    midsig = round(length(sig) / 2);
    idx = round(midsig - 2 * fs * 60 : midsig + 2 * fs * 60);
    idx2 = iis.peakPos > idx(1) & iis.peakPos < idx(end);
    plot(tstamps(idx) / 60, sig(idx), 'k')
    axis tight
    hold on
    scatter(iis.peakPos(idx2) / fs / 60,...
    iis.peakPower(idx2), '*');
    ylabel('Voltage [mV]')
    xlabel('Time [m]')
    xticks(round([midsig / fs / 60 - 2, midsig / fs / 60 + 2]))
    set(gca, 'TickLength', [0 0])
    box off
    title('IIS')
    
    % iis waveforms
    subplot(3, 4, 3)
    plot(wvstamps * 1000, seg')
    ylabel('Voltage [mV]')
    xlabel('Time [ms]')
    axis tight
    xticks([-marg, 0, marg] * 1000);
    set(gca, 'TickLength', [0 0])
    box off
    title('Spike waveform')
       
    % mean + std waveform
    axes('Position',[.542 .71 .09 .07])
    box on
    stdshade(seg, 0.5, 'k', wvstamps)
    axis tight
    set(gca, 'TickLength', [0 0], 'YTickLabel', [], 'XTickLabel', [],...
        'XColor', 'none', 'YColor', 'none', 'Color', 'none')
    title(sprintf('n = %d', size(seg, 1)));
    box off
    
    % iis cwt
    subplot(3, 4, 4)
    imagesc(wvstamps, f, abs(cfs))
    axis xy
    set(gca, 'YScale', 'log')
    % caxis([0.0001 1])
    origSize = get(gca, 'Position');
    colorbar
    set(gca, 'Position', origSize);
    hold on
    plot(wvstamps, coi, 'w', 'LineWidth', 2)
    ylim([min(f) max(f)]);
    xlim([wvstamps(1) wvstamps(end)]);
    yticks([ceil(min(f)), 50, 100, max(f)])
    xticks([min(wvstamps), 0, max(wvstamps)])
    ylabel('Frequency [Hz]')
    xlabel('Time [ms]')
    set(gca, 'TickLength', [0 0])
    box off
    title('Average Scalogram')
    
    % weidth histogram
    if ~isempty(iis.maxf)
        subplot(3, 4, 7)
        h = histogram(log10(iis.spkw), 20, 'Normalization', 'Probability');
        h.EdgeColor = 'none';
        h.FaceColor = 'k';
        h.FaceAlpha = 1;
        xlabel('Spike Width [ms]')
        ylabel('Probability [%]')
        set(gca, 'TickLength', [0 0])
        box off
        title('IIS width')
    end
    
    % amplitude histogram
    subplot(3, 4, 8)
    h = histogram(log10(iis.peakPower), 30, 'Normalization', 'Probability');
    h.EdgeColor = 'none';
    h.FaceColor = 'k';
    h.FaceAlpha = 1;
    xlabel('Peak voltage [log(uV)]')
    ylabel('Probability [%]')
    set(gca, 'TickLength', [0 0])
    box off
    title('Amplitude Distribution')
    
    % max frequency and amplitude vs. time
    subplot(3, 4, 7)
    scatter(iis.peakPos / fs / 60, iis.peakPower, 2, 'b', 'filled');
    axis tight
    x = xlim;
    l2 = lsline;
    set(l2, 'color', 'b')
    l2.LineWidth = 3;
    xlim(x);
    ylabel('Amplitude [mV]')
    if ~isempty(iis.maxf)
        yyaxis right
        scatter(iis.peakPos / fs / 60, iis.maxf, 2, 'k', 'filled');
        l1 = lsline;
        ylabel('Frequency [Hz]')
        set(l1, 'color', 'k')
        l1.LineWidth = 3;
    end
    axis tight
    xlabel('Time [m]')
    set(gca, 'TickLength', [0 0])
    box off
    title('Amplitude. vs. Time')
    
    % CCH 5s
    subplot(3, 4, 11)
    plotCCG('ccg', ccg, 't', tccg / 1000, 'basepath', basepath,...
        'saveFig', false, 'c', {'k'});
    xlabel('Time [s]')
    title('Autocorrelogram')
    
    % CCH 10s
    subplot(3, 4, 12)
    plotCCG('ccg', ccg2, 't', tccg2 / 1000, 'basepath', basepath,...
        'saveFig', false, 'c', {'k'});
    xlabel('Time [s]')
    title('Autocorrelogram')
    
    if saveFig
        figname = [basename '_IIS'];
        export_fig(figname, '-tif', '-transparent')
        % savePdf(figname, basepath, ff)
    end
end
end

% EOF

