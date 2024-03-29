function swv = spkwv_metrics(varargin)

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
%   loadSpk     logical. load spikes from .spk files
% 
% OUTPUT
%   swv         struct
%      
% DEPENDENCIES
% 
% TO DO LIST
%   # add tail slope (Torrado Pacheco et al., Neuron, 2021). TP threshold
%   in that article was ~0.4 ms (done)
%   # snip across all tetrode channels to present in cell explorer
%   # handle inverted spikes, possible by flipping before calculating
%   params
%
% 08 apr 19 LH      updates: 
% 14 may 20 LH      added upsampling by fft
% 12 dec 21 LH      snip raw from dat
% 29 dec 21 LH      added time to repolarization
% 10 feb 23 LH      handle cases were minimum is at end of wv

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
addOptional(p, 'loadSpk', false, @islogical);

parse(p, varargin{:})
basepath    = p.Results.basepath;
wv          = p.Results.wv;
fs          = p.Results.fs;
saveVar     = p.Results.saveVar;
forceA      = p.Results.forceA;
loadSpk     = p.Results.loadSpk;

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
datfile = fullfile(basepath, [basename, '.dat']);

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

% size of one data point in bytes
precision = 'int16';
nbytes = class2bytes(precision);

% number of spikes to snip per cluster
spks2snip = 10000;

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
% that the spikes struct of CE is still necassary to obtain the maximum
% amplitude channel before resnipping.

if isempty(wv)
        
    if exist(swvRawFile, 'file')
        load(swvRawFile)
        
    elseif ~loadSpk
        % memory map to datfile
        info = dir(datfile);
        nsamps = info.bytes / nbytes / nchans;
        m = memmapfile(datfile, 'Format', {precision, [nchans, nsamps] 'mapped'});
        raw = m.Data;
        
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
        
        [swv_raw, ~] = snipFromBinary('stamps', spktimes, 'fname', datfile,...
            'win', win, 'nchans', nchans, 'ch', ch, 'align_peak', 'min',...
            'precision', precision, 'rmv_trend', 6, 'saveVar', false,...
            'l2norm', false, 'raw', raw);
        swv_raw = cellfun(@squeeze, swv_raw, 'UniformOutput', false);
        
        if saveVar
            save(swvRawFile, 'swv_raw')
        end
        
        clear raw
        clear m
    else
        
        load(spkFile);
        grps = unique(spikes.shankID);
        cnt = 1;
        
        % load from spk file
        for igrp = grps
            spk = loadNS('datatype', 'spk', 'session', session, 'grpid', igrp);
            clu = loadNS('datatype', 'clu', 'session', session, 'grpid', igrp);
            uclu = unique(clu);
            for iclu = 1 : length(uclu)
                if uclu(iclu) ~= 0 && uclu(iclu) ~= 1
                    
                    % randomly select a subset of spikes 
                    cluidx = find(clu == uclu(iclu));
                    randidx = randperm(length(cluidx), min([length(cluidx), spks2snip]));
                    cluidx = cluidx(randidx)
                    
                    % get channel of maximum amplitude
                    tmp = spk(:, :, cluidx);
                    [~, maxCh] = max(range(mean(tmp, 3), 2));
                    swv_raw{cnt} = squeeze(tmp(maxCh, :, :));
                    cnt = cnt + 1;
                end
            end
        end

        if saveVar
            save(swvRawFile, 'swv_raw')
        end
    end
    
    % l2norm
    for iunit = 1 : length(swv_raw)
        swv_raw{iunit} = swv_raw{iunit} ./ vecnorm(swv_raw{iunit}, 2, 1);
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
imin = nan(1, nunits);

for iunit = 1 : nunits
    
    % upsampled waveform
    w = wv_interp(iunit, :);    
        
    % spike width by inverse of max frequency in spectrum (stark et al.,
    % 2013)
    [cfs, f, ~] = cwt(w, 'FilterBank', fb);
    [~, ifreq] = max(abs(squeeze(cfs)), [], 1);
    maxf = f(ifreq(round(length(w) / 2)));
    spkw(iunit) = 1000 / maxf;

    % general waveform params
    [minVal, imin(iunit)] = min(w);                             % trough
    [maxVal_post, imax_post] = max(w(imin(iunit) + 1 : end));   % peak after trough
    imax_post = imax_post + imin(iunit);
    [maxVal_pre, ~] = max(w(1 : imin(iunit) - 1));              % peak before trough
    [minVal_post, imin_post] = min(w(imax_post + 1 : end));     % min after peak
    if isempty(minVal_post)
        minVal_post = w(end);
        imin_post = length(w);
    end
    imin_post = imin_post + imax_post;
    
    if isempty(maxVal_post)
        continue
    end

    % amplitudes
    ampTp(iunit) = maxVal_post - minVal;                        % amplitude trough to after peak
    ampTail(iunit) = maxVal_post - minVal_post;                 % amplitude tail
    
    % trough-to-peak time (artho et al., 2004) and asymmetry (Sirota et
    % al., 2008)
    if ~isempty(imax_post)
        tp(iunit) = (imax_post - imin(iunit)) * 1000 / nfs;      % samples to ms
        if ~isempty(maxVal_pre)
            asym(iunit) = (maxVal_post - maxVal_pre) / (maxVal_post + maxVal_pre);
        end
    else
        warning('waveform may be corrupted')
    end
    
    % slope peak to end (Torrado Pacheco et al., Neuron, 2021)
    slopeTail(iunit) = ampTail(iunit) / (imin_post - imax_post);
    
    % slope trough to peak (no reference)
    slopeTp(iunit) = ampTp(iunit) / (imax_post - imin(iunit));
    
    % half peak width (Medrihan et al., 2017)
    wu = w / maxVal_post;
    th1 = find(wu(imin(iunit) : imax_post) < 0.5);
    if any(th1)
        th1 = imax_post - th1(end);
        th2 = find(wu(imax_post + imin(iunit) : end) < 0.5);
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
    
    % complex spike index (McHugh 1996), defined as percentage of first lag
    % isi that fall between 3 ms and 15 ms and whose second spike is
    % smaller than the first. could not find a code reference so
    % implemented manually. csi requires knowing the amp of each spike.
    % here I snip ~20k spikes which is ok so long as they are contineous in
    % time and not randomely seleceted.
    
    
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
swv.imin = imin;

if saveVar       
    save(swvFile, 'swv')
    
    % cell metrics
    if exist(cmFile, 'file')
        load(cmFile)
%         cell_metrics.waveforms.wv = num2cell(swv.wv, 2)';
        cell_metrics.waveforms.wv = num2cell(wv_interp, 2)';
%         cell_metrics.waveforms.wv_std = num2cell(swv.wv_std, 2)';
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
