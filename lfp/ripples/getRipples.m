function ripp = getRipples(varargin)

% Detects ripples from lfp: based in part on bz_FindRipples, Liu et al.,
% Nat. Comm., 2022, and TL
% 
% INPUT:
%   basepath        path to recording {pwd}
%   sig             numeric vec of lfp data. if empty will be loaded from
%                   basename.lfp according to rippCh.  
%   emg             numeric of emg data. must be the same sampling
%                   frequency as sig. can be acc.mag. used for
%                   normalization and exclusion of artifacts. must be
%                   compatible with recWin. if not will try to fix this by
%                   assuming emg represents the entire recording
%   fs          	numeric. sampling frequency of lfp file / data. 
%                   if empty will be extracted from session info (ce
%                   format)
%   rippCh          numeric vec. channels to load and average from lfp
%                   file {1}. 
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
%   binary2epochs
%   lfpFilter
%   Sync (buzcode)
%   SyncMap (buzcode)
%   PlotColorMap (buzcode)
%
% TO DO LIST:
%   # finish graphics (done)
%   # stats (done)
%   # rate (done)
%   # exclusion by emg noise (done)
%   # exlude active periods when normalizing signal (done)
%   # allow user to input sig directly instead of loading from binary (done)
%   # exclusion by spiking activity (irrelevant - can be done manually)
%   # improve routine to select best ripple channel automatically (done)
%   # add subroutine to determine bandpass, threshold, and duration limits
%   automatically by performing a quick first detection
%   # add graphics for visualizing detection independently of neuroscope
%   # implement batch proccessing for memory
%   # use multiple channels for detection
%   # use phase difference between tetrodes to exclude artifacts
%   # add option to exclude ripples detected on a "bad" channel 
%
% 02 dec 21 LH     updates:
% 11 jan 23        implementad remarks from Liu et al., Nat. Comm., 2022
%                  separated getRippleSpks and graphics    
% 25 jan 23        removed subroutine to detect best ripples channel and
%                  normalize signal by nrem

% BATCH PROCESSING: a 24 hr flat binary (int16) file of 19 channels at 1250
% Hz is ~4 GB so with high end computers, batch processing is not critical.
% detection includes the signal, filtered signal, amplitude (from hilbert),
% and maybe also phase (from hilbert), and the smooth amplitude, all as double. 

% detect ripples on each channel separately (good and
% bad). requires filtering, amplitude, smooth amplitude. goes until
% binary2epochs. than if overlap keep (good) or remove (bad).
% next, average the start / end times from all channels. 
% next, average the signal and recalculate, amp, phase, etc.
% start exclusion. 

% getRipples_exclusion: 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'sig', [], @isnumeric);
addOptional(p, 'emg', [], @isnumeric);
addOptional(p, 'recWin', [0 Inf]);
addOptional(p, 'rippCh', [], @isnumeric);
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
fs          = p.Results.fs;
graphics    = p.Results.graphics;
saveVar     = p.Results.saveVar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% files
cd(basepath)
[~, basename] = fileparts(basepath);
ssfile = fullfile(basepath, [basename '.sleep_states.mat']);
sessionfile = fullfile(basepath, [basename, '.session.mat']);
rippfile = fullfile(basepath, [basename, '.ripp.mat']);
clufile = fullfile(basepath, [basename, '.ripp.clu.1']);
resfile = fullfile(basepath, [basename, '.ripp.res.1']);

% load session info
if ~exist(sessionfile, 'file')
    session = CE_sessionTemplate(pwd, 'viaGUI', false,...
        'forceDef', false, 'forceL', false, 'saveVar', false);
else
    load(sessionfile)
end

% session params
spkgrp = session.extracellular.spikeGroups.channels;
spkch = sort([spkgrp{:}]);
nchans = session.extracellular.nChannels;
fsSpks = session.extracellular.sr;
if isempty(fs)
    fs = session.extracellular.srLfp;
end

% threshold of stds above sig_amp. initial detection, peak power,
% contineous power, artifact power
thr = [1.5, 2.5, 2, 200];
% thr = [0.5, 1, 1.5, 200];
% thr = [0.5, 1, 0.5, 200];

% detection params. limDir refers to min, max, and inter-ripple
% duration. 4th element refers to the amount of time the power must be
% above thr(3). all values in [ms]
limDur = [20, 300, 20, 10]; 
limDur = round(limDur / 1000 * fs);
passband = [100 300];
binsizeRate = 60;           % binsize for calculating ripple rate [s]
emgThr = 75;                % exclude ripples that occur when emg > thr

fprintf('\ngetting ripples for %s\n', basename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(sig)    

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
% prepare signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('preparing signals...\n')

% prepare emg. assumes emg and sig are of the same sample frequency. if the
% lengths are incompatible, assumes emg represents the entire recording and
% clips it according to recWin
if length(emg) > length(sig)  
    sigIdx = recWin(1) * fs + 1 : recWin(2) * fs;
    emg = emg(sigIdx);
end

% filter lfp data in ripple band
sig_filt = filterLFP(sig, 'fs', fs, 'type', 'butter', 'dataOnly', true,...
    'order', 5, 'passband', passband, 'graphics', false);

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

% clear memory
clear sig_unwrapped h 
        
% -------------------------------------------------------------------------
% detection signal
dtctMet = 1;            % 1 = TL; 2 = BZ
switch dtctMet
    case 1              % TL detection (smoothed amplitude)
        % smooth amplitude. note TL averages the smoothed amplitude from a few
        % channels.
        sig_dtct = smoothdata(sig_amp, 'gaussian', round(5 / (1000 / 1250)));       
        sig_dtct = (sig_dtct - mean(sig_dtct)) / std(sig_dtct);
        
    case 2              % BZ detection (normalized squared signal)
        sig_dtct = sig_filt .^ 2;
        winFilt = ones(11, 1) / 11;
        shift = (length(winFilt) - 1) / 2;
        [sig_dtct, z] = filter(winFilt, 1, sig_dtct);
        sig_dtct = [sig_dtct(shift + 1 : end, :); z(1 : shift, :)];
        sig_dtct = (sig_dtct - mean(sig_dtct)) / std(sig_dtct);
            
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find ripples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find epochs with power > low threshold. correct for durations
epochs = binary2epochs('vec', sig_dtct > thr(1), 'minDur', limDur(1),...
    'maxDur', limDur(2), 'interDur', limDur(3));
nepochs = size(epochs, 1);

% discard ripples that occur during high emg 
if ~isempty(emg)
    
    % calc emg rms throughout the recording in 4s bins
    emg_binLen = fs * 4;
    emg_nbins = floor(length(emg) / emg_binLen);
    emg_rms = log(rms(reshape(emg(1 : emg_nbins * emg_binLen), emg_binLen, emg_nbins)));

    % calc emg rms per epoch
    for iepoch = 1 : nepochs
        emgRipp(iepoch) = log(rms(emg(epochs(iepoch, 1) : epochs(iepoch, 2))));
    end
    discard_idx = emgRipp > prctile(emg_rms, emgThr);

%     find threshold from the bimodal distribution of emg
%         [~, cents] = kmeans(emg_rms(:), 2);
%         emgThr = mean(cents);
%         discard_idx = emgRipp > emgThr;
    epochs(discard_idx, :) = [];
    emgRipp(discard_idx) = [];
    nepochs = size(epochs, 1);
    fprintf('After emg exclusion: %d events\n', nepochs)
end

% discard ripples with a peak power < high threshold
clear discard_idx
peakPowNorm = zeros(size(epochs, 1), 1);
for iepoch = 1 : size(epochs, 1)
    peakPowNorm(iepoch) = max(sig_dtct(epochs(iepoch, 1) : epochs(iepoch, 2)));
end
discard_idx = peakPowNorm < thr(2);
epochs(discard_idx, :) = [];
peakPowNorm(discard_idx) = [];
nepochs = size(epochs, 1);
fprintf('After peak power: %d events\n', nepochs)

% discard ripples with a peak power > atrtifact threshold
clear discard_idx
discard_idx = peakPowNorm > thr(4);
epochs(discard_idx, :) = [];
peakPowNorm(discard_idx) = [];
nepochs = size(epochs, 1);
fprintf('After artifact power: %d events\n', nepochs)

% discard ripples that do not maintain peak power for the min duration
discard_idx = false(1, nepochs);
for iepoch = 1 : size(epochs, 1)
    aboveThr = sig_dtct(epochs(iepoch, 1) : epochs(iepoch, 2)) > thr(3);
    powPnts = strfind([0 aboveThr'], [0 ones(1, limDur(4))]);
    if isempty(powPnts)
        discard_idx(iepoch) = true;
    end
end
epochs(discard_idx, :) = [];
peakPowNorm(discard_idx) = [];
nepochs = size(epochs, 1);
fprintf('After contineous power: %d events\n', nepochs)

% clear memory
clear sig_dtct

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ripp stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% change to maxPow and minPow

% find negative peak position for each ripple
peakPos = zeros(size(epochs, 1), 1);
peakPow = zeros(size(epochs, 1), 1);
for iepoch = 1 : size(epochs, 1)
    [peakPow(iepoch), peakPos(iepoch)] =...
        min(sig_filt(epochs(iepoch, 1) : epochs(iepoch, 2)));
    peakPos(iepoch) = peakPos(iepoch) + epochs(iepoch, 1) - 1;
end

% convert idx to seconds
peakPos = peakPos / fs;
epochs = epochs / fs;

% -------------------------------------------------------------------------
% maps
ripp.maps.durWin = [-75 75] / 1000;
nbinsMap = floor(fs * diff(ripp.maps.durWin) / 2) * 2 + 1; % must be odd
centerBin = ceil(nbinsMap / 2);
[r, i] = Sync([tstamps sig], peakPos, 'durations', ripp.maps.durWin);
ripp.maps.raw = SyncMap(r, i, 'durations', ripp.maps.durWin,...
    'nbins', nbinsMap, 'smooth', 0);
[r, i] = Sync([tstamps sig_filt], peakPos, 'durations', ripp.maps.durWin);
ripp.maps.ripp = SyncMap(r, i, 'durations', ripp.maps.durWin,...
    'nbins', nbinsMap, 'smooth', 0);
[f, i] = Sync([tstamps, sig_freq], peakPos, 'durations', ripp.maps.durWin);
ripp.maps.freq = SyncMap(f, i, 'durations', ripp.maps.durWin,...
    'nbins', nbinsMap, 'smooth', 0);
[a, i] = Sync([tstamps sig_phase], peakPos, 'durations', ripp.maps.durWin);
ripp.maps.phase = SyncMap(a, i,'durations', ripp.maps.durWin,...
    'nbins', nbinsMap, 'smooth', 0);
[p, i] = Sync([tstamps sig_amp], peakPos, 'durations', ripp.maps.durWin);
ripp.maps.amp = SyncMap(p, i,'durations', ripp.maps.durWin,...
    'nbins', nbinsMap, 'smooth', 0);

% clear memory
clear sig_freq sig_amp sig_phase tstamps sig_filt

% -------------------------------------------------------------------------
% more stats
ripp.maxFreq = max(ripp.maps.freq, [], 2);
ripp.peakFreq = ripp.maps.freq(:, centerBin);
ripp.peakAmp = ripp.maps.amp(:, centerBin);
ripp.dur = epochs(:, 2) - epochs(:, 1);

% acg and correlations
[ripp.acg.data, ripp.acg.t] = CCG(peakPos,...
    ones(length(peakPos), 1), 'binSize', 0.01);
ripp.corr.amp_freq = corrcoef(ripp.peakAmp, ripp.peakFreq);
ripp.corr.dur_freq = corrcoef(ripp.dur, ripp.peakFreq);
ripp.corr.dur_amp = corrcoef(ripp.dur, ripp.peakAmp);

% rate of ripples
epochs             = epochs + recWin(1);  
peakPos            = peakPos + recWin(1);
[ripp.rate.rate, ripp.rate.binedges, ripp.rate.tstamps] =...
    times2rate(peakPos, 'binsize', binsizeRate, 'winCalc', recWin,...
    'c2r', true);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finalize, save, and graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ripp.info.rippCh        = rippCh;
ripp.info.limDur        = limDur;
ripp.info.recWin        = recWin;
ripp.info.runtime       = datetime(now, 'ConvertFrom', 'datenum');
ripp.info.thr           = thr;
ripp.info.binsizeRate   = binsizeRate;
ripp.info.fs            = fs;
ripp.info.passband      = passband;
ripp.epochs             = epochs;  
ripp.peakPos            = peakPos;
ripp.peakPow            = peakPow;
ripp.peakPowNorm        = peakPowNorm;
   
if saveVar      
    save(rippfile, 'ripp')

    % create ns files for visualization with neuroscope. note currently
    % this only works if recWin(1) = 0
    nepochs = size(ripp.epochs, 1);

    res = round([ripp.epochs(:, 1); ripp.epochs(:, 2); ripp.peakPos] * fsSpks);
    [res, sort_idx] = sort(res);
    fid = fopen(resfile, 'w');
    fprintf(fid, '%d\n', res);
    rc = fclose(fid);
    if rc == -1
        error('failed to close res file')
    end

    clu = [ones(nepochs, 1); ones(nepochs, 1) * 2; ones(nepochs, 1) * 3];
    clu = clu(sort_idx);
    fid = fopen(clufile, 'w');
    fprintf(fid, '%d\n', 3);
    fprintf(fid, '%d\n', clu);
    rc = fclose(fid);
    if rc == -1
        error('failed to close clu file')
    end
end

if graphics 
    plot_ripples(ripp, 'basepath', basepath, 'saveFig', true)
end
 

end

% EOF