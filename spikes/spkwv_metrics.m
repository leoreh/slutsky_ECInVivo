function swv = spkwv_metrics(varargin)
% SPKWV_METRICS Calculates various waveform parameters for cell classification.
%
% SUMMARY:
% This function calculates waveform parameters typically used to separate
% regular spiking (RS) and fast spiking (FS) cells. It can work with either
% pre-computed waveforms or extract them from raw data files. The function
% supports both .dat and .spk file formats, automatically choosing the
% appropriate source based on file availability.
%
% METHODOLOGY:
% The function first loads or extracts waveforms using spkwv_load. Then it
% calculates various metrics including:
%   - Trough-to-peak duration and amplitude
%   - Spike width (inverse of max frequency)
%   - Waveform asymmetry
%   - Half peak width
%   - Repolarization time
%   - Tail slope and amplitude
%   - Peak duration and ratio metrics
%
% Waveforms are upsampled (default 20x) using either spline interpolation
% or FFT for improved temporal resolution of the metrics.
%
% INPUT (Optional Key-Value Pairs):
%   basepath    (char) Full path to recording folder. If empty, uses current
%               directory. {pwd}
%   wv          (numeric matrix) Pre-computed waveforms matrix of units
%               (rows) x samples (columns). If empty, waveforms will be
%               extracted from raw data files. {[]}
%   fs          (numeric) Sampling frequency in Hz. Used to determine
%               waveform length (typically 1.6 ms). If empty, will be loaded
%               from session.mat. {[]}
%   flgSave     (logical) Flag to save the output struct to a .swv_metrics.mat
%               file. {true}
%   flgForce    (logical) Flag to force analysis even if output file exists.
%               {false}
%
% OUTPUT:
%   swv         (struct) Structure containing waveform metrics:
%               .tp          - Trough-to-peak duration (ms)
%               .tpAmp       - Trough-to-peak amplitude
%               .tpRatio     - Peak-to-trough ratio
%               .tpSlope     - Slope from trough to peak
%               .spkw        - Spike width (ms)
%               .asym        - Waveform asymmetry
%               .hpk         - Half peak width (ms)
%               .rtau        - Repolarization time (ms)
%               .tailSlope   - Slope from peak to end
%               .tailAmp     - Amplitude of tail
%               .peakDuration - Duration between peaks (ms)
%               .oneMinusLeftPeak - 1 - left peak value
%
% DEPENDENCIES:
%   spkwv_load, cwtfilterbank
%
% HISTORY:
%   08 apr 19 LH      Initial version
%   14 may 20 LH      Added upsampling by fft
%   12 dec 21 LH      Added dat file support
%   29 dec 21 LH      Added time to repolarization
%   10 feb 23 LH      Handle cases where minimum is at end of wv
%   Aug 2024          Updated documentation and streamlined code
%
% TO DO:
% complex spike index (McHugh 1996), defined as percentage of first lag
% isi that fall between 3 ms and 15 ms and whose second spike is
% smaller than the first. could not find a code reference so
% implemented manually. csi requires knowing the amp of each spike.
% here I snip ~20k spikes which is ok so long as they are contineous in
% time and not randomely seleceted.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'wv', []);
addOptional(p, 'fs', [], @isnumeric);
addOptional(p, 'flgSave', true, @islogical);
addOptional(p, 'flgForce', false, @islogical);

parse(p, varargin{:})
basepath    = p.Results.basepath;
wv          = p.Results.wv;
fs          = p.Results.fs;
flgSave     = p.Results.flgSave;
flgForce    = p.Results.flgForce;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file names
[~, basename] = fileparts(basepath);
swvFile = fullfile(basepath, [basename, '.swv_metrics.mat']);
sessionFile = fullfile(basepath, [basename, '.session.mat']);

% check if already analyzed
if exist(swvFile, 'file') && ~flgForce
    load(swvFile)
    return
end

% get params from session info
if exist(sessionFile, 'file')
    load(sessionFile);
    if isempty(fs)
        fs = session.extracellular.sr;
    end
end

% waveform params
spklength = ceil(1.6 * 10^-3 * fs);  % spike wave is 1.6 ms
upsamp = 20;            % upsample factor for waveforms
upsamp_met = 'spline';  % method for upsampling. can also be 'fft'.

% create wavelet filter
nfs = fs * upsamp;
fb = cwtfilterbank('SignalLength', spklength * upsamp, 'VoicesPerOctave', 32,...
    'SamplingFrequency', nfs, 'FrequencyLimits', [1 fs / 4]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(wv)
    
    % load raw waveforms using spkwv_load
    swv_raw = spkwv_load('basepath', basepath, 'fs', fs);

    if ~isempty(swv_raw)
        % perform l2norm using MATLAB's vecnorm
        swv_raw = cellfun(@(x) x ./ vecnorm(x, 2, 1),...
            swv_raw, 'UniformOutput', false);

        % get mean and std
        wv = cellfun(@(x) [mean(x, 2)]', swv_raw, 'UniformOutput', false);
        wv = cell2mat(wv');
        wv_std = cellfun(@(x) [std(x, [], 2)]', swv_raw, 'UniformOutput', false);
        wv_std = cell2mat(wv_std');
    
    else
        % backup: load from spikes struct
        spkFile = fullfile(basepath, [basename, '.spikes.cellinfo.mat']);
        load(spkFile, 'spikes');
        swv_raw = cellfun(@(x) x ./ vecnorm(x', 2, 1),...
            spikes.rawWaveform', 'UniformOutput', false);

        wv = cell2mat(swv_raw);
        wv_std = cell2mat(spikes.rawWaveform_std');
    end
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

% initialize
nunits = size(wv, 1);
tp = nan(1, nunits);
spkw = nan(1, nunits);
hpk = nan(1, nunits);
asym = nan(1, nunits);
rtau = nan(1, nunits);
tpSlope = zeros(1, nunits);
tailSlope = zeros(1, nunits);
tpAmp = nan(1, nunits);
tailAmp = nan(1, nunits);
imin = nan(1, nunits);
tpRatio = nan(1, nunits);
peakDuration = nan(1, nunits);
oneMinusLeftPeak = nan(1, nunits);
imax_pre = nan(1, nunits);
inverted = false(1, nunits);

for iunit = 1 : nunits

    % upsampled waveform
    w = wv_interp(iunit, :);

    % spike width by inverse of max frequency in spectrum (stark et al.,
    % 2013)
    [cfs, f, ~] = cwt(w, 'FilterBank', fb);
    [~, ifreq] = max(abs(squeeze(cfs)), [], 1);
    maxf = f(ifreq(round(length(w) / 2)));
    spkw(iunit) = 1000 / maxf;

    % check if waveforn inverted
    if abs(max(w)) > abs(min(w))
        inverted(iunit) = true;
    end

    % general waveform params
    [minVal, imin(iunit)] = min(w);                                     % trough
    [maxVal_post, imax_post] = max(w(imin(iunit) + 1 : end));           % peak after trough
    imax_post = imax_post + imin(iunit);
    [maxVal_pre, imax_pre(iunit)] = max([w(1 : imin(iunit) - 1), 1]);   % peak before trough
    [minVal_post, imin_post] = min(w(imax_post + 1 : end));             % min after peak
    if isempty(minVal_post)
        minVal_post = w(end);
        imin_post = length(w);
    end
    imin_post = imin_post + imax_post;

    if isempty(maxVal_post)
        continue
    end

    % amplitudes
    tpAmp(iunit) = maxVal_post - minVal;                        % amplitude trough to after peak
    tailAmp(iunit) = maxVal_post - minVal_post;                 % amplitude tail

    % trough-to-peak time (artho et al., 2004) and asymmetry (Sirota et
    % al., 2008)
    if ~isempty(imax_post)
        tp(iunit) = (imax_post - imin(iunit)) * 1000 / nfs;      % samples to ms
        if ~isempty(maxVal_pre)
            asym(iunit) = (maxVal_post - maxVal_pre) / (maxVal_post + maxVal_pre);

            % Calculate new metrics
            tpRatio(iunit) = abs(maxVal_post) / abs(minVal);  % peak-to-trough ratio
            peakDuration(iunit) = (imax_post - imax_pre(iunit)) * 1000 / nfs;  % duration between peaks in ms
            oneMinusLeftPeak(iunit) = 1 - maxVal_pre;  % 1 - left peak value
        end
    else
        warning('waveform may be corrupted')
    end

    % slope peak to end (Torrado Pacheco et al., Neuron, 2021)
    tailSlope(iunit) = tailAmp(iunit) / (imin_post - imax_post);

    % slope trough to peak (no reference)
    tpSlope(iunit) = tpAmp(iunit) / (imax_post - imin(iunit));

    % half peak width (Medrihan et al., 2017).
    [posPk, ~, posWdth] = findpeaks(w, nfs, 'SortStr', 'descend');
    [negPk, negLocs, negWdth] = findpeaks(-w, nfs, 'SortStr', 'descend');
    pkWdth = [posWdth, negWdth];
    [~, pkIdx] = max([posPk, negPk]);
    hpk(iunit) = pkWdth(pkIdx) * 1000;

    % time for repolarization (Ardid et al., J. Neurosci., 2015;
    % https://github.com/LofNaDI). this fails for most of our
    % cells.
    decayVal = maxVal_post - 0.25 * tpAmp(iunit);
    rtau_idx = nearest(w(imax_post + 1 : end), decayVal);
    if ~isempty(rtau_idx)
        rtau(iunit) = x_time(imax_post + rtau_idx) - x_time(imax_post);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

swv.info.runtime = datetime("now");
swv.info.upsamp_met = 'spline';
swv.info.upsamp = upsamp;
swv.wv = wv;
swv.wv_std = wv_std;
swv.tp = tp;
swv.tpAmp = tpAmp;
swv.tpRatio = tpRatio;
swv.tpSlope = tpSlope;
swv.spkw = spkw;
swv.asym = asym;
swv.hpk = hpk;
swv.rtau = rtau;
swv.tailSlope = tailSlope;
swv.tailAmp = tailAmp;
swv.imin = imin;
swv.peakDuration = peakDuration;
swv.oneMinusLeftPeak = oneMinusLeftPeak;
swv.inverted = inverted;

if flgSave
    save(swvFile, 'swv')
end

end     % EOF
