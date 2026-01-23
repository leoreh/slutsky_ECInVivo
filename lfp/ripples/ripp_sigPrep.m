function rippSig = ripp_sigPrep(lfp, fs, varargin)
% RIPP_SIGPREP Prepares LFP signals for ripple detection.
%
%   rippSig = RIPP_SIGPREP(lfp, fs, varargin)
%
%   SUMMARY:
%       Processes raw LFP to extract features needed for detection:
%       1.  Bandpass filtering (default 100-300Hz).
%       2.  Hilbert transform for Amplitude and Frequency.
%       3.  Generation of a "Detection Signal" (e.g., Smoothed Squared Energy).
%       4.  Z-Scoring of the detection signal (Adaptive or NREM-based).
%
%   INPUTS:
%       lfp         - (Vec) Raw LFP signal (Voltages).
%       fs          - (Num) Sampling frequency [Hz].
%       varargin    - Parameter/Value pairs:
%           'detectMet' - (Num) Method ID for detection signal (default: 3).
%                           1: Smoothed Amplitude (Hilbert Env).
%                           2: Smoothed Squared Amplitude.
%                           3: Smoothed TEO (Teager Energy Operator).
%                           4: Rectified + Lowpass (4th order Butter).
%           'passband'  - (Vec) Filtering range [min max] (Hz). Default: [100 300].
%           'zMet'      - (Char) Normalization method:
%                           'adaptive' : Moving average/std (10s window).
%                           'nrem'     : Global mean/std derived from NREM epochs.
%           'nremTimes' - (Mat) [N x 2] NREM start/end times (Required for 'nrem').
%
%   OUTPUTS:
%       rippSig     - (Struct) Processed signals:
%           .lfp        - Original Raw LFP.
%           .filt       - Bandpass filtered LFP.
%           .amp        - Amplitude envelope (Hilbert).
%           .freq       - Instantaneous frequency (Hilbert).
%           .z          - Final Z-scored detection signal.
%
%   DEPENDENCIES:
%       filterLFP
%
%   HISTORY:
%       Updated: 23 Jan 2026

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
% 1. Filter LFP
rippSig.lfp = lfp;
rippSig.filt = filterLFP(lfp, 'fs', fs, 'type', 'butter', 'dataOnly', true,...
    'order', 5, 'passband', passband, 'graphics', false);

% 2. Analytic Signal (Hilbert Transform)
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

% 4. Generate Base Detection Signal
switch detectMet
    case 1 % Smoothed Amplitude
        baseSignal = rippSig.amp;
        winSmooth = round(0.005 * fs);

    case 2 % Smoothed Power (Squared)
        baseSignal = rippSig.filt .^ 2;
        winSmooth = round(0.010 * fs);

    case 3 % Smoothed Teager Energy Operator (TEO)
        sPad = [rippSig.filt(1); rippSig.filt; rippSig.filt(end)];
        baseSignal = sPad(2:end-1).^2 - sPad(1:end-2) .* sPad(3:end);
        baseSignal(baseSignal < 0) = 0; % TEO should be positive
        winSmooth = round(0.005 * fs);

    case 4 % Rectified + Lowpass Filter
        rectSig = abs(rippSig.filt);
        lpFreq = 55;
        [b_lp, a_lp] = butter(4, lpFreq / (fs / 2), 'low');
        baseSignal = filtfilt(b_lp, a_lp, rectSig);
        winSmooth = 0; % No additional smoothing
end

if winSmooth > 0
    baseSignal = smoothdata(baseSignal, 'gaussian', winSmooth);
end

% 5. Z-Score Normalization
switch zMet
    case 'adaptive'
        % Local Moving Average/Std
        movLen = round(10 * fs);
        mu = movmean(baseSignal, movLen);
        sigma = movstd(baseSignal, movLen);

        % Prevent division by zero/noise
        stdFloor = 1e-6 * mean(sigma, 'omitnan');
        sigma(sigma < stdFloor) = stdFloor;

    case 'nrem'
        % Global Mean/Std from NREM epochs only

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
