function rippSig = ripp_sigPrep(lfp, fs, varargin)
% RIPP_SIGPREP Prepares LFP signals for ripple detection.
%
%   rippSig = RIPP_SIGPREP(lfp, fs, varargin)
%
%   SUMMARY:
%       Processes raw LFP to extract features needed for detection:
%       1.  Bandpass filtering.
%       2.  Hilbert transform for Amplitude and Frequency.
%       3.  Generation of a "Detection Signal" (e.g., Smoothed Squared Energy).
%       4.  Z-Scoring of the detection signal (Adaptive or NREM-based).
%
%   INPUTS:
%       lfp         - (Vec) Raw LFP signal (Voltages).
%       fs          - (Num) Sampling frequency [Hz].
%       varargin    - Parameter/Value pairs:
%           'detectMet' - (Num) Method ID for detection signal.
%                           1: Smoothed Amplitude (Hilbert Env).
%                           2: Smoothed Squared Amplitude.
%                           3: Smoothed TEO (Teager Energy Operator).
%                           4: Rectified + Lowpass (4th order Butter).
%           'passband'  - (Vec) Filtering range [min max] (Hz).
%           'zMet'      - (Char) Normalization method:
%                           'adaptive' : Moving average/std (10s window).
%                           'nrem'     : Global mean/std from NREM epochs.
%                           'nremBg'   : NREM mean/std with candidate ripples
%                                        removed (signal-independent baseline;
%                                        keeps a real group rate difference from
%                                        being normalized away). Falls back to
%                                        'adaptive' when NREM is absent.
%           'nremTimes' - (Mat) [N x 2] NREM start/end times (for 'nrem'/'nremBg').
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
        winSmooth = round(0.015 * fs);

    case 4 % Rectified, clipped, & Lowpass Filter
        rectSig = abs(rippSig.filt);
        
        % Calculate "Robust" SD based on MAD for clipping
        clipVal = 4 * median(rectSig) / 0.6745;
        rectSig(rectSig > clipVal) = clipVal;

        % Smooth the clipped signal
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
        % Global Mean/Std from NREM epochs only. Falls back to 'adaptive'
        % when no NREM samples are available (missing / empty sleep states),
        % so detection still runs on sessions without vigilance scoring.
        mask = false(size(baseSignal));
        if ~isempty(nremTimes)
            nremSamp = round(nremTimes * fs) + 1;
            nSamples = length(baseSignal);
            for iBout = 1:size(nremSamp, 1)
                idxStart = max(1, nremSamp(iBout, 1));
                idxEnd = min(nSamples, nremSamp(iBout, 2));
                if idxStart <= idxEnd
                    mask(idxStart:idxEnd) = true;
                end
            end
        end

        if any(mask)
            mu = mean(baseSignal(mask), 'omitnan');
            sigma = std(baseSignal(mask), 'omitnan');
        else
            warning('ripp_sigPrep:noNrem', ...
                ['No NREM samples for ''nrem'' z-scoring; ', ...
                'falling back to ''adaptive''.']);
            movLen = round(10 * fs);
            mu = movmean(baseSignal, movLen);
            sigma = movstd(baseSignal, movLen);
            stdFloor = 1e-6 * mean(sigma, 'omitnan');
            sigma(sigma < stdFloor) = stdFloor;
        end

    case 'nremBg'
        % Signal-INDEPENDENT baseline: NREM background with the candidate
        % ripples removed. The plain 'nrem' baseline includes the ripples, so
        % its SD scales with ripple power - which inflates the threshold in a
        % high-ripple animal and can normalize away a genuine group difference
        % in ripple rate. Referencing the noise floor (ripples excluded) makes
        % the threshold depend on the background only, so a real group effect
        % survives. Use this to test whether per-animal z-scoring is hiding an
        % effect (see the detection review, P3). Falls back like 'nrem'.
        mask = false(size(baseSignal));
        nSamples = length(baseSignal);
        if ~isempty(nremTimes)
            nremSamp = round(nremTimes * fs) + 1;
            for iBout = 1:size(nremSamp, 1)
                idxStart = max(1, nremSamp(iBout, 1));
                idxEnd = min(nSamples, nremSamp(iBout, 2));
                if idxStart <= idxEnd
                    mask(idxStart:idxEnd) = true;
                end
            end
        end

        if any(mask)
            mu0 = mean(baseSignal(mask), 'omitnan');
            sd0 = std(baseSignal(mask), 'omitnan');
            if sd0 == 0, sd0 = 1; end

            % Drop supra-threshold samples (candidate ripples), dilated +/-25 ms,
            % then re-estimate the baseline on what remains.
            isEvt  = (baseSignal - mu0) / sd0 > 2;
            isEvt  = movmax(double(isEvt), round(0.050 * fs)) > 0;
            bgMask = mask & ~isEvt;
            if nnz(bgMask) < 0.1 * nnz(mask)    % too little left -> keep all NREM
                bgMask = mask;
            end
            mu    = mean(baseSignal(bgMask), 'omitnan');
            sigma = std(baseSignal(bgMask), 'omitnan');
            if sigma == 0, sigma = 1; end
        else
            warning('ripp_sigPrep:noNrem', ...
                ['No NREM samples for ''nremBg'' z-scoring; ', ...
                'falling back to ''adaptive''.']);
            movLen = round(10 * fs);
            mu = movmean(baseSignal, movLen);
            sigma = movstd(baseSignal, movLen);
            stdFloor = 1e-6 * mean(sigma, 'omitnan');
            sigma(sigma < stdFloor) = stdFloor;
        end

    otherwise
        error('Unknown zMet: %s', zMet);
end

rippSig.z = (baseSignal - mu) ./ sigma;

end         % EOF
