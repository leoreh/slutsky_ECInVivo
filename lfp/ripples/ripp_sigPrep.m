function rippSig = ripp_sigPrep(lfp, fs, varargin)
% RIPP_SIGPREP Prepares signals for ripple detection and analysis.
%
% INPUT:
%   lfp             - (Vec) Raw LFP signal.
%   fs              - (Num) Sampling rate in Hz.
%   varargin        - 'detectAlt' (default 3), 'passband' ([100 300])
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
%   3. Calculate detection signal (based on detectAlt) and Z-score it.

% Parameters
p = inputParser;
addRequired(p, 'lfp', @isnumeric);
addRequired(p, 'fs', @isnumeric);
addParameter(p, 'detectAlt', 3, @isnumeric);
addParameter(p, 'passband', [100 300], @isnumeric);
parse(p, lfp, fs, varargin{:});

detectAlt = p.Results.detectAlt;
passband = p.Results.passband;

% Filter LFP for detection
sigFilt = filterLFP(lfp, 'fs', fs, 'type', 'butter', 'dataOnly', true,...
    'order', 5, 'passband', passband, 'graphics', false);

% Hilbert Analysis
sigH = hilbert(sigFilt);
sigAmp = abs(sigH);
sigPhase = angle(sigH);
sigUnwrapped = unwrap(sigPhase);

% Instantaneous Frequency
timestamps = (0:length(lfp)-1)' / fs;
dt = 1/fs;
d0 = diff(medfilt1(sigUnwrapped, 12)) ./ dt;
t_diff = timestamps(1:end-1) + dt/2;
d1 = interp1(t_diff, d0, timestamps(2:end-1), 'linear', 'extrap');
sigFreq = [d0(1); d1; d0(end)] / (2 * pi);

% Detection Signal Prep
switch detectAlt
    case 1 % Smoothed Amp
        baseSignal = sigAmp;
        winSmooth = round(0.005 * fs);
    case 2 % Smoothed Sq
        baseSignal = sigFilt .^ 2;
        winSmooth = round(0.010 * fs);
    case 3 % Smoothed TEO
        sPad = [sigFilt(1); sigFilt; sigFilt(end)];
        baseSignal = sPad(2:end-1).^2 - sPad(1:end-2) .* sPad(3:end);
        baseSignal(baseSignal < 0) = 0;
        winSmooth = round(0.005 * fs);
    case 4 % Rectified + Lowpass
        rectSig = abs(sigFilt);
        lpFreq = 55;
        [b_lp, a_lp] = butter(4, lpFreq / (fs / 2), 'low');
        baseSignal = filtfilt(b_lp, a_lp, rectSig);
        winSmooth = 0;
end

if winSmooth > 0
    baseSignal = smoothdata(baseSignal, 'gaussian', winSmooth);
end

% Z-Score (Adaptive Threshold)
movLen = round(10 * fs);
localMean = movmean(baseSignal, movLen);
localStd = movstd(baseSignal, movLen);
stdFloor = 1e-6 * mean(localStd, 'omitnan');
localStd(localStd < stdFloor) = stdFloor;
z = (baseSignal - localMean) ./ localStd;

% Output Struct
rippSig = struct();
rippSig.lfp = lfp;
rippSig.filt = sigFilt;
rippSig.amp = sigAmp;
rippSig.freq = sigFreq;
rippSig.z = z;
rippSig.sigPhase = sigPhase;

end         % EOF
