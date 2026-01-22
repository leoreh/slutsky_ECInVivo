function rippSig = ripp_sigPrep(lfp, fs, varargin)
% RIPP_SIGPREP Prepares signals for ripple detection and analysis.
%
% INPUT:
%   lfp             - (Vec) Raw LFP signal.
%   fs              - (Num) Sampling rate in Hz.
%   varargin        - 'detectMet' (default 3), 'passband' ([100 300])
%
% OUTPUT:
%   sig         - Structure with fields:
%       .lfp        - Raw LFP
%       .filt       - Filtered LFP
%       .amp        - Amplitude envelope
%       .freq       - Instantaneous frequency
%       .z          - Z-scored detection signal
%
% LOGIC:
%   1. Filter LFP (Butterworth 100-300Hz by default)
%   2. Calculate envelope/frequency via Hilbert
%   3. Calculate detection signal (based on detectMet) and Z-score it.

% Parameters
p = inputParser;
addRequired(p, 'lfp', @isnumeric);
addRequired(p, 'fs', @isnumeric);
addParameter(p, 'detectMet', 3, @isnumeric);
addParameter(p, 'passband', [100 300], @isnumeric);
addParameter(p, 'zMet', 'adaptive', @ischar);
addParameter(p, 'nremTimes', [], @isnumeric);
parse(p, lfp, fs, varargin{:});

detectMet = p.Results.detectMet;
passband = p.Results.passband;
zMet = p.Results.zMet;
nremTimes = p.Results.nremTimes;

% Filter LFP for detection
% Filter LFP for detection
rippSig.lfp = lfp;
rippSig.filt = filterLFP(lfp, 'fs', fs, 'type', 'butter', 'dataOnly', true,...
    'order', 5, 'passband', passband, 'graphics', false);

% Hilbert Analysis
sigH = hilbert(rippSig.filt);
rippSig.amp = abs(sigH);
rippSig.sigPhase = angle(sigH);
sigUnwrapped = unwrap(rippSig.sigPhase);

% Instantaneous Frequency
timestamps = (0:length(lfp)-1)' / fs;
dt = 1/fs;
d0 = diff(medfilt1(sigUnwrapped, 12)) ./ dt;
t_diff = timestamps(1:end-1) + dt/2;
d1 = interp1(t_diff, d0, timestamps(2:end-1), 'linear', 'extrap');
rippSig.freq = [d0(1); d1; d0(end)] / (2 * pi);

% Detection Signal Prep
switch detectMet
    case 1 % Smoothed Amp
        baseSignal = rippSig.amp;
        winSmooth = round(0.005 * fs);
    case 2 % Smoothed Sq
        baseSignal = rippSig.filt .^ 2;
        winSmooth = round(0.010 * fs);
    case 3 % Smoothed TEO
        sPad = [rippSig.filt(1); rippSig.filt; rippSig.filt(end)];
        baseSignal = sPad(2:end-1).^2 - sPad(1:end-2) .* sPad(3:end);
        baseSignal(baseSignal < 0) = 0;
        winSmooth = round(0.005 * fs);
    case 4 % Rectified + Lowpass
        rectSig = abs(rippSig.filt);
        lpFreq = 55;
        [b_lp, a_lp] = butter(4, lpFreq / (fs / 2), 'low');
        baseSignal = filtfilt(b_lp, a_lp, rectSig);
        winSmooth = 0;
end

if winSmooth > 0
    baseSignal = smoothdata(baseSignal, 'gaussian', winSmooth);
end

% Z-Score
switch zMet
    case 'adaptive'
        % Moving average/std
        movLen = round(10 * fs);
        mu = movmean(baseSignal, movLen);
        sigma = movstd(baseSignal, movLen);
        stdFloor = 1e-6 * mean(sigma, 'omitnan');
        sigma(sigma < stdFloor) = stdFloor;

    case 'nrem'
        % Global mean/std based on NREM times

        mask = false(size(baseSignal));
        nremSamp = round(nremTimes * fs) + 1;
        nSamples = length(baseSignal);

        for iBout = 1:size(nremSamp, 1)
            idxStart = max(1, nremSamp(iBout, 1));
            idxEnd = min(nSamples, nremSamp(iBout, 2));
            if idxStart <= idxEnd
                mask(idxStart:idxEnd) = true;
            end
        end

        mu = mean(baseSignal(mask), 'omitnan');
        sigma = std(baseSignal(mask), 'omitnan');

    otherwise
        error('Unknown zMet: %s', zMet);
end

rippSig.z = (baseSignal - mu) ./ sigma;

end         % EOF
