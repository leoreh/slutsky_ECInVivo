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
%   in that article was ~0.4 ms
%
% 08 apr 19 LH      updates: 
% 14 may 20 LH      added upsampling by fft
% 12 dec 21 LH      snip raw from dat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
swvRawFile = fullfile(basepath, [basename, '.swv.mat']);
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
spks2snip = 10000;

% waveform params
if ~isempty(wv)
    spklength = size(wv, 2);
else
    spklength = ceil(1.6 * 10^-3 * fs);  % spike wave is 1.6 ms
end
win = [-(floor(spklength / 2) - 1) floor(spklength / 2)];   
upsamp = 10;            % upsample factor for waveforms
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
% if mean waveforms were not given, resnip from the raw dat file a
% random fraction of spikes for each cluster. this is much superior to the
% extraction done by CE (accurate detrending instead of
% filtering). note that the spikes struct of CE is currently still
% necassary for the maximum amplitude channel

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
        [wv_all, ~] = snipFromBinary('stamps', spktimes, 'fname', datname,...
            'win', win, 'nchans', nchans, 'ch', ch, 'align_peak', 'min',...
            'precision', 'int16', 'rmv_trend', 6, 'saveVar', false,...
            'l2norm', false);
        wv_all = cellfun(@squeeze, wv_all, 'UniformOutput', false);
        
        if saveVar
            save(swvRawFile, 'wv_all')
        end
        
    end
    
    % get mean and std
    wv = cellfun(@(x) [mean(x, 2)]', wv_all, 'UniformOutput', false);
    wv = cell2mat(wv');
    wv_std = cellfun(@(x) [std(x, [], 2)]', wv_all, 'UniformOutput', false);
    wv_std = cell2mat(wv_std');

end

% interpolate
x_orig = linspace(0, 1, size(wv, 2));
x_upsamp = linspace(0, 1, size(wv, 2) * upsamp);
switch upsamp_met
    case 'spline'
        wv_interp = [interp1(x_orig, wv', x_upsamp, 'spline', nan)]';
    case 'fft'
        wv_interp = interpft(wv, upsamp * size(wv, 2), 2); % upsample in the frequency domain
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\ncalculating waveform parameters\n\n')
tic

% initialize
nunits = size(wv, 1);
tp = nan(1, nunits);
spkw = nan(1, nunits);
hpk = nan(1, nunits);
asym = nan(1, nunits);

for iunit = 1 : nunits
    
    w = wv_interp(iunit, :);
    
    % ---------------------------------------------------------------------
    % trough-to-peak time (artho et al., 2004) and asymmetry (Sirota et
    % al., 2008)
    [~, minpos] = min(w);
    [maxval, ~] = max(w(1 : minpos - 1));   
    [maxvalpost, maxpost] = max(w(minpos + 1 : end));               
    if ~isempty(maxpost)
        tp(iunit) = maxpost;
        if ~isempty(maxval)
            asym(iunit) = (maxvalpost - maxval) / (maxvalpost + maxval);
        end
    else
        warning('waveform may be corrupted')
        tp(iunit) = NaN;
        asym(iunit) = NaN;
    end
    
    % ---------------------------------------------------------------------
    % half peak width (Medrihan et al., 2017)
    wu = w / maxvalpost;
    th1 = find(wu(minpos : maxpost + minpos) < 0.5);
    if any(th1)
        th1 = maxpost - th1(end);
        th2 = find(wu(maxpost + minpos : end) < 0.5);        
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
    
    % add Ardid et al., J. Neurosci., 2015 (https://github.com/LofNaDI)
    
    %  spike width by inverse of max frequency in spectrum (stark et al., 2013)
    [cfs, f, ~] = cwt(w, 'FilterBank', fb);
    [~, ifreq] = max(abs(squeeze(cfs)), [], 1);
    maxf = f(ifreq(round(length(w) / 2)));
    spkw(iunit) = 1000 / maxf;
end

% samples to ms
tp = tp * 1000 / nfs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

swv.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
swv.info.upsamp_met = 'spline'; 
swv.info.upsamp = upsamp; 
swv.wv = wv;
swv.wv_std = wv_std;
swv.tp = tp;
swv.spkw = spkw;
swv.asym = asym;
swv.hpk = hpk;

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
        save(cmFile, 'cell_metrics')
    end
end

fprintf('\nthat took %.1f minutes\n', toc / 60)

return

% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paths used to compare CE results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepaths{1} = 'G:\RA\hDLX_Gq_WT2\200820_bslDay1';
basepaths{2} = 'G:\RA\hDLX_Gq_Tg\210820_bslDay2Raw2';
basepaths{3} = 'D:\Data\lh86\lh86_210301_072600';
basepaths{4} = 'G:\lh81\lh81_210207_045300';
cell_metrics = CellExplorer('basepaths', basepaths);

