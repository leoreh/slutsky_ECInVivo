function ripp = getRipples(varargin)

% based on bz_FindRipples. several detection params were changed based on
% manual visualization of results. 


% INPUT:
%   basepath        path to recording {pwd}
%   fs          	numeric. sampling frequency of lfp file. if empty will
%                   be extracted from session info (ce format)
%   ch              numeric vec. channels to load and average from lfp
%                   file {1}.
%   recWin          numeric 2 element vector. time to anlayse in recording
%                   [s]. start of recording is marked by 0. {[0 Inf]}
%   graphics        logical. plot graphics {true} or not (false)
%   saveVar         logical. save variables (update spikes and save su)
%
% OUTPUT:
%   ripp            struct
%
% DEPENDENCIES:
%   binary2epochs
%   lfpFilter
%   Sync (buzcode)
%   SyncMap (buzcode)
%   PlotColorMap (buzcode)
%
% TO DO LIST:
%   finish graphics
%   stats (done)
%   rate
%
% 02 dec 21 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'recWin', [0 Inf]);
addOptional(p, 'ch', 1, @isnumeric);
addOptional(p, 'fs', [], @isnumeric);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
basepath    = p.Results.basepath;
recWin      = p.Results.recWin;
ch          = p.Results.ch;
fs          = p.Results.fs;
graphics    = p.Results.graphics;
saveVar     = p.Results.saveVar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = pwd;
win = ones(11, 1) / 11;     % for moving average
thr = [1 2.5];              % threshold of stds above the sig
limDur = [20, 150, 30];     % min, max, and inter dur limits for ripples [ms]
passband = [100 200];

% load session info
[~, basename] = fileparts(basepath);
sessionName = [basename, '.session.mat'];
if ~exist(sessionName, 'file')
    session = CE_sessionTemplate(pwd, 'viaGUI', false,...
        'force', true, 'saveVar', true);
else
    load(sessionName)
end

nchans = session.extracellular.nChannels;
if isempty(fs)
    fs = session.extracellular.srLfp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data
sig = double(bz_LoadBinary([basename, '.lfp'], 'duration', diff(recWin),...
    'frequency', fs, 'nchannels', nchans, 'start', recWin(1),...
    'channels', ch, 'downsample', 1));
sig = mean(sig, 2);

% filter
sig_filt = filterLFP(sig, 'fs', fs, 'type', 'butter', 'dataOnly', true,...
    'order', 3, 'passband', passband, 'graphics', false);

% square
nss = sig_filt.^2;

% move mean and correct shift
[nss, z] = filter(win, 1, nss);
shift = (length(win) - 1) / 2;
nss = [nss(shift + 1 : end, :); z(1 : shift, :)];

% standardize
nss = (nss - mean(nss)) / std(nss); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find ripples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find epochs with power > low threshold. correct for durations
limDur = limDur / 1000 * fs;
epochs = binary2epochs('vec', nss > thr(1), 'minDur', limDur(1),...
    'maxDur', limDur(2), 'interDur', limDur(3));

% discard ripples with a peak power < high threshold
peakPowNorm = zeros(size(epochs, 1), 1);
for iepoch = 1 : size(epochs, 1)
    peakPowNorm(iepoch) = max(nss(epochs(iepoch, 1) : epochs(iepoch, 2)));
end
dicard_idx = peakPowNorm < thr(2);
epochs(dicard_idx, :) = [];
peakPowNorm(dicard_idx) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find negative peak position for each ripple
peakPos = zeros(size(epochs, 1), 1);
peakPow = zeros(size(epochs, 1), 1);
for iepoch = 1 : size(epochs, 1)
    [peakPow(iepoch), peakPos(iepoch)] =...
        min(sig(epochs(iepoch, 1) : epochs(iepoch, 2)));
    peakPos(iepoch) = peakPos(iepoch) + epochs(iepoch, 1) - 1;
end

% convert idx to seconds
nepochs = size(epochs, 1);
peakPos = peakPos / fs;
epochs = epochs / fs;

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

ripp.peakFreq = ripp.maps.freq(:, centerBin);
ripp.peakAmp = ripp.maps.amp(:,centerBin);
ripp.dur = epochs(:, 2) - epochs(:, 1);

% acg and correlations
[ripp.acg.data, ripp.acg.t] = CCG(peakPos,...
    ones(length(peakPos), 1), 'binSize', 0.01);
ripp.corr.amp_freq = corrcoef(ripp.peakAmp, ripp.peakFreq);
ripp.corr.dur_freq = corrcoef(ripp.dur, ripp.peakFreq);
ripp.corr.dur_amp = corrcoef(ripp.dur, ripp.peakAmp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    fh = figure;
    durPlot = [-50 50] / 1000;
    x = durPlot(1) : diff(durPlot) / nepochs : durPlot(2);
    histBins = 200;
    
    sb1 = subplot(3, 3, 1);
    ripp_idx = randperm(nepochs, 1000);
    plot(((1:nbins)' - ceil(nbins / 2)) / nbins * diff(durPlot), ripp.maps.ripp(ripp_idx, :)', 'k');
    xlabel('Time [s]')
    
    sb2 = subplot(3, 3, [2, 3]);
    
    sb4 = subplot(3, 3, 4);
    PlotColorMap(ripp.maps.freq, 1, 'bar','on', 'cutoffs', [100 250], 'x', x);
    xlabel('Ripple Frequency');
    
    sb5 = subplot(3, 3, 5);
    PlotColorMap(ripp.maps.amp, 1, 'bar','on', 'x', x);
    xlabel('Ripple Amplitude');
    
    sb6 = subplot(3, 3, 6);
    plotCCG(ripp.acg.data, ripp.acg.t);
    xlabel('Time [ms]')
    ylabel('Rate')
    
    sb7 = subplot(3, 3, 7);
    h = histogram(ripp.peakFreq, histBins, 'Normalization', 'probability');
    h.FaceColor = 'k';
    h.EdgeColor = 'none';
    xlabel('Peak Frequency [Hz]')
    ylabel('Probability')
    
    sb8 = subplot(3, 3, 8);
    h = histogram(ripp.peakAmp, histBins, 'Normalization', 'probability');
    h.FaceColor = 'k';
    h.EdgeColor = 'none';
    xlabel('Peak Amp')
    ylabel('Probability')
    
    sb9 = subplot(3, 3, 9);
    h = histogram(ripp.dur * 1000, histBins, 'Normalization', 'probability');
    h.FaceColor = 'k';
    h.EdgeColor = 'none';
    xlabel('Ripple Duration [ms]')
    ylabel('Probability')
    
    figpath = fullfile(basepath, 'graphics');
    mkdir(figpath)
    figname = fullfile(figpath, sprintf('%s_ripples', basename));
    export_fig(figname, '-tif', '-transparent', '-r300')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finalize and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ripp.info.ch = ch;
ripp.info.limDur = limDur;
ripp.info.recWin = recWin;
ripp.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
ripp.info.thr = thr;
ripp.epochs         = epochs + recWin(1);  % convert to seconds
ripp.peakPos        = peakPos;
ripp.peakPow        = peakPow;
ripp.peakPowNorm    = peakPowNorm;

if saveVar      
    save([basepath, filesep, basename, '.ripp.mat'], 'ripp')
    
    % create ns files for visualization with neuroscope
    fs_dat = session.extracellular.sr;
    fs_dat = 24414.0625;
    
    res = [ripp.epochs(:, 1); ripp.epochs(:, 2); ripp.peakPos] * fs_dat;
    [res, sort_idx] = sort(round(res));
    fid = fopen([basename, '.ripp.res.1'], 'w');
    fprintf(fid, '%d\n', res);
    rc = fclose(fid);
   
    clu = [ones(nepochs, 1); ones(nepochs, 1) * 2; ones(nepochs, 1) * 3];
    clu = clu(sort_idx);
    fid = fopen([basename, '.ripp.clu.1'], 'w');
    fprintf(fid, '%d\n', 3);
    fprintf(fid, '%d\n', clu);
    rc = fclose(fid);
  
end

end

% EOF