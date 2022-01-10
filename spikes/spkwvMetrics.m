function swv = spkwvMetrics(varargin)

% calculates various waveform parameters typically used to separate RS and
% FS cells
% 
% INPUT
%   basepath    recording session path {pwd}
%   wv          matrix of units (rows) x samples (columns). if empty,
%               waveforms will be extracted from the raw dat file. for
%               example: wv = cat(1, spikes.rawWaveform{spikes.su})'
%   fs          sampling frequency. used only to determine the waveform
%               length (which is typically 1.6 ms)
%   saveVar     save variable {1}.
%   forceA      logical. force analysis even if struct file exists {false}
% 
% OUTPUT
%   swv         struct
%      
% DEPENDENCIES
% 
% TO DO LIST
%   # add tail slope (Torrado Pacheco et al., Neuron, 2021). TP threshold
%   in that article was ~0.4 ms (done)
%
% 08 apr 19 LH      updates: 
% 14 may 20 LH      added upsampling by fft
% 12 dec 21 LH      snip raw from dat
% 29 dec 21 LH      added time to repolarization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'wv', []);
addOptional(p, 'fs', [], @isnumeric);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'forceA', false, @islogical);

parse(p, varargin{:})
basepath    = p.Results.basepath;
wv          = p.Results.wv;
fs          = p.Results.fs;
saveVar     = p.Results.saveVar;
forceA      = p.Results.forceA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file names
[~, basename] = fileparts(basepath);
cmFile = fullfile(basepath, [basename, '.cell_metrics.cellinfo.mat']);
swvFile = fullfile(basepath, [basename, '.swv_metrics.mat']);
swvRawFile = fullfile(basepath, [basename, '.swv_raw.mat']);
sessionFile = fullfile(basepath, [basename, '.session.mat']);
spkFile = fullfile(basepath, [basename, '.spikes.cellinfo.mat']);

% check if already analyzed 
if exist(swvFile, 'file') && ~forceA
    load(swvFile)
    return
end

% get params from session info
if exist(sessionFile, 'file')
    load(sessionFile);
    if isempty(fs)
        fs = session.extracellular.sr;
    end
    nchans = session.extracellular.nChannels;
end

% number of spikes to snip per cluster
spks2snip = 20000;

% waveform params
if ~isempty(wv)
    spklength = size(wv, 2);
else
    spklength = ceil(1.6 * 10^-3 * fs);  % spike wave is 1.6 ms
end
win = [-(floor(spklength / 2) - 1) floor(spklength / 2)];   
upsamp = 20;            % upsample factor for waveforms
upsamp_met = 'spline';  % method for upsampling. can also be 'fft'.

% create wavelet filter
nfs = fs * upsamp;
fb = cwtfilterbank('SignalLength', spklength * upsamp, 'VoicesPerOctave', 32,...
    'SamplingFrequency', nfs, 'FrequencyLimits', [1 fs / 4]);

% initialize
wv_std = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if mean waveforms were not given, resnip from the raw dat file a random
% number of spikes for each cluster. this is much superior to the
% extraction done by CE (accurate detrending instead of filtering). note
% that the spikes struct of CE is still necassary to obtain for the maximum
% amplitude channel before resnipping.

if isempty(wv)
        
    if exist(swvRawFile, 'file')
        load(swvRawFile)
        
    else

        % load spikes and session info struct
        load(spkFile);
        ch = num2cell(spikes.maxWaveformCh1);
        
        % select a random fraction of spikes
        fn = @(x) randperm(length(x), min([length(x), spks2snip]));
        spkidx = cellfun(fn, spikes.ts, 'UniformOutput', false);
        
        % clip spktimes
        for iunit = 1 : length(spikes.ts)
            spktimes{iunit} = spikes.ts{iunit}(spkidx{iunit});
        end
        
        datname = fullfile(basepath, [basename, '.dat']);
        [swv_raw, ~] = snipFromBinary('stamps', spktimes, 'fname', datname,...
            'win', win, 'nchans', nchans, 'ch', ch, 'align_peak', 'min',...
            'precision', 'int16', 'rmv_trend', 6, 'saveVar', false,...
            'l2norm', false);
        swv_raw = cellfun(@squeeze, swv_raw, 'UniformOutput', false);
        
        if saveVar
            save(swvRawFile, 'swv_raw')
        end
        
    end
    
    % get mean and std
    wv = cellfun(@(x) [mean(x, 2)]', swv_raw, 'UniformOutput', false);
    wv = cell2mat(wv');
    wv_std = cellfun(@(x) [std(x, [], 2)]', swv_raw, 'UniformOutput', false);
    wv_std = cell2mat(wv_std');

end

% interpolate
x_orig = linspace(0, 1, size(wv, 2));
x_upsamp = linspace(0, 1, size(wv, 2) * upsamp);
x_time = [1 : size(wv, 2) * upsamp] / nfs * 1000;
switch upsamp_met
    case 'spline'
        wv_interp = [interp1(x_orig, wv', x_upsamp, 'spline', nan)]';
    case 'fft'
        wv_interp = interpft(wv, upsamp * size(wv, 2), 2); % upsample in the frequency domain
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\ncalculating waveform parameters\n')

% initialize
nunits = size(wv, 1);
tp = nan(1, nunits);
spkw = nan(1, nunits);
hpk = nan(1, nunits);
asym = nan(1, nunits);
rtau = nan(1, nunits);
slopeTp = zeros(1, nunits);
slopeTail = zeros(1, nunits);
ampTp = nan(1, nunits);
ampTail = nan(1, nunits);

for iunit = 1 : nunits
    
    % general waveform params
    w = wv_interp(iunit, :);
    [minVal, imin] = min(w);                                % trough
    [maxVal_post, imax_post] = max(w(imin + 1 : end));      % peak after trough
    imax_post = imax_post + imin;
    [maxVal_pre, ~] = max(w(1 : imin - 1));                 % peak before trough
    [minVal_post, imin_post] = min(w(imax_post + 1 : end)); % min after peak
    if isempty(minVal_post)
        minVal_post = w(end);
        imin_post = length(w);
    end
    imin_post = imin_post + imax_post;
    ampTp(iunit) = maxVal_post - minVal;                    % amplitude trough to after peak
    ampTail(iunit) = maxVal_post - minVal_post;             % amplitude tail
    
    % trough-to-peak time (artho et al., 2004) and asymmetry (Sirota et
    % al., 2008)
    if ~isempty(imax_post)
        tp(iunit) = imax_post * 1000 / nfs;      % samples to ms
        if ~isempty(maxVal_pre)
            asym(iunit) = (maxVal_post - maxVal_pre) / (maxVal_post + maxVal_pre);
        end
    else
        warning('waveform may be corrupted')
    end
    
    % slope peak to end (Torrado Pacheco et al., Neuron, 2021)
    slopeTail(iunit) = ampTail(iunit) / (imin_post - imax_post);
    
    % slope trough to peak (no reference)
    slopeTp(iunit) = ampTp(iunit) / (imax_post - imin);

    
    % half peak width (Medrihan et al., 2017)
    wu = w / maxVal_post;
    th1 = find(wu(imin : imax_post) < 0.5);
    if any(th1)
        th1 = imax_post - th1(end);
        th2 = find(wu(imax_post + imin : end) < 0.5);
        if any(th2)
            th2 = th2(1);
        else
            th2 = th1;
        end
    else
        th1 = NaN;
        th2 = NaN;
    end
    hpk(iunit) = (th1 + th2) * 1000 / nfs;
    
    % time for repolarization (Ardid et al., J. Neurosci., 2015;
    % https://github.com/LofNaDI). this fails for most of our
    % cells.
    decayVal = maxVal_post - 0.25 * ampTp(iunit);
    rtau_idx = nearest(w(imax_post + 1 : end), decayVal);
    if ~isempty(rtau_idx)
        rtau(iunit) = x_time(imax_post + rtau_idx) - x_time(imax_post);
    end
    
    % spike width by inverse of max frequency in spectrum (stark et al., 2013)
    [cfs, f, ~] = cwt(w, 'FilterBank', fb);
    [~, ifreq] = max(abs(squeeze(cfs)), [], 1);
    maxf = f(ifreq(round(length(w) / 2)));
    spkw(iunit) = 1000 / maxf;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

swv.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
swv.info.upsamp_met = 'spline'; 
swv.info.upsamp = upsamp; 
swv.info.spks2snip = spks2snip;
swv.wv = wv;
swv.wv_std = wv_std;
swv.tp = tp;
swv.spkw = spkw;
swv.asym = asym;
swv.hpk = hpk;
swv.rtau = rtau;
swv.slopeTp = slopeTp;
swv.slopeTail = slopeTail;
swv.ampTp = ampTp;
swv.ampTail = ampTail;

if saveVar       
    save(swvFile, 'swv')
    
    % cell metrics
    if exist(cmFile, 'file')
        load(cmFile)
        cell_metrics.waveforms.wv = num2cell(swv.wv, 2)';
        cell_metrics.waveforms.filt = num2cell(swv.wv, 2)';
        cell_metrics.waveforms.filt_std = num2cell(swv.wv_std, 2)';
        cell_metrics.swv_spkw = swv.spkw;
        cell_metrics.swv_tp = swv.tp;
        cell_metrics.swv_asym = swv.asym;
        cell_metrics.swv_hpk = swv.hpk;
        cell_metrics.swv_rtau = swv.rtau;
        cell_metrics.swv_slopeTp = swv.slopeTp;
        cell_metrics.swv_slopeTail = swv.slopeTail;
        cell_metrics.swv_ampTp = swv.ampTp;
        cell_metrics.swv_ampTail = swv.ampTail;
        save(cmFile, 'cell_metrics')
    end
end

fprintf('\nthat took %.1f minutes\n', toc / 60)


return

% EOF
