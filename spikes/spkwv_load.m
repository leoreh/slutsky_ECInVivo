function swv_raw = spkwv_load(varargin)
% SPKWV_LOAD Extracts raw waveforms from dat or spk files for spike sorting and analysis.
%
% SUMMARY:
% This function extracts raw spike waveforms from either .dat or .spk files.
% It can load waveforms directly from a pre-saved .swv_raw.mat file, or
% extract them from the raw data files. For .dat files, it uses memory
% mapping for efficient data access and includes detrending. For .spk files,
% it loads waveforms directly from the spike files. The function supports
% both single-channel and multi-channel recordings.
%
% METHODOLOGY:
% The function first checks for an existing .swv_raw.mat file. If not found,
% it either:
%   1. Extracts waveforms from .dat file:
%      - Memory maps the .dat file for efficient access
%      - Uses spike times from .spikes.cellinfo.mat
%      - Snips waveforms around each spike using snipFromBinary
%      - Applies detrending and L2 normalization
%   2. Loads waveforms from .spk files:
%      - Loads spike data per shank
%      - Selects channel with maximum amplitude
%      - Applies L2 normalization
%
% INPUT (Optional Key-Value Pairs):
%   basepath    (char) Full path to recording folder. If empty, uses current
%               directory. {pwd}
%   fs          (numeric) Sampling frequency in Hz. Used to determine
%               waveform length (typically 1.6 ms). If empty, will be loaded
%               from session.mat. {[]}
%   flgSave     (logical) Flag to save the raw waveforms to a .swv_raw.mat
%               file. {true}
%
% OUTPUT:
%   swv_raw     (cell array) Cell array of raw waveforms per unit. Each cell
%               contains a matrix of spikes (rows) x samples (columns).
%               Waveforms are L2 normalized.
%
% DEPENDENCIES:
%   snipFromBinary, loadNS, class2bytes
%
% HISTORY:
%   12 dec 21 LH      Initial version with dat file support
%   29 dec 21 LH      Added spk file support
%   Aug 2024          Updated documentation and added flgSave option

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'fs', [], @isnumeric);
addOptional(p, 'flgSave', true, @islogical);

parse(p, varargin{:})
basepath    = p.Results.basepath;
fs          = p.Results.fs;
flgSave     = p.Results.flgSave;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file names
[~, basename] = fileparts(basepath);
sessionFile = fullfile(basepath, [basename, '.session.mat']);
spkFile = fullfile(basepath, [basename, '.spikes.cellinfo.mat']);
datfile = fullfile(basepath, [basename, '.dat']);
swvRawFile = fullfile(basepath, [basename, '.swv_raw.mat']);

% check if raw waveforms already extracted 
if exist(swvRawFile, 'file')
    load(swvRawFile)
    return
end

% check if dat file exists, if not force loadSpk
if ~exist(datfile, 'file')
    loadSpk = true;
    warning('dat file not found, loading from spk files instead')
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
spklength = ceil(1.6 * 10^-3 * fs);  % spike wave is 1.6 ms
win = [-(floor(spklength / 2) - 1) floor(spklength / 2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if mean waveforms were not given, resnip from the raw dat file a random
% number of spikes for each cluster. this is much superior to the
% extraction done by CE (accurate detrending instead of filtering). note
% that the spikes struct of CE is still necassary to obtain the maximum
% amplitude channel before resnipping.

if ~loadSpk
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

    if flgSave
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
                cluidx = cluidx(randidx);

                % get channel of maximum amplitude
                tmp = spk(:, :, cluidx);
                [~, maxCh] = max(range(mean(tmp, 3), 2));
                swv_raw{cnt} = squeeze(tmp(maxCh, :, :));
                cnt = cnt + 1;
            end
        end
    end

    if flgSave
        save(swvRawFile, 'swv_raw')
    end
end

end % EOF
