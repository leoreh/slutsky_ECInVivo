function ripp = ripp_params(rippSig, ripp)
% RIPP_PARAMS Calculates physiological parameters for ripple events.
%
%   ripp = RIPP_PARAMS(rippSig, ripp)
%
%   SUMMARY:
%       Computes "Best Practice" parameters for each event using the
%       pre-processed signals.
%       1. Duration: Exact time difference (End - Start).
%       2. Peak Amplitude: Envelope amplitude at the peak index.
%       3. Mean Frequency: Average instantaneous frequency (Hz) in fixed window.
%       4. Event Frequency: Average freq over detection duration (Hz).
%       5. Total Energy: Sum of squared filtered signal (Integrated Power).
%       6. Skewness: Center of Mass of the envelope relative to peak (in ms).
%       7. Peak Frequency: Spectral peak of the whitened event PSD (Hz). The
%          Hilbert instantaneous frequency (.freq/.freqEvent) is biased LOW by
%          the 1/f slope; .freqPeak removes the aperiodic tilt and reports the
%          true oscillation frequency. Additive - the Hilbert measures are kept.
%
%   INPUTS:
%       rippSig     - (Struct) Signal structure (requires .lfp, .filt, .amp, .freq).
%       ripp        - (Struct) Event structure with .times and .peakTime.
%                              Must contain .info.fs.
%
%   OUTPUTS:
%       ripp        - (Struct) Updated structure with new fields:
%                       .amp       (N x 1) [uV]
%                       .freq      (N x 1) [Hz] (Fixed window, Hilbert)
%                       .freqEvent (N x 1) [Hz] (Full Duration, Hilbert)
%                       .freqPeak  (N x 1) [Hz] (Whitened spectral peak)
%                       .peakProm  (N x 1) [ratio] (Whitened peak height above
%                                  the 1/f background; ~1 for a broadband
%                                  transient, >>1 for a true ripple)
%                       .energy    (N x 1) [uV^2]
%                       .dur       (N x 1) [ms]
%                       .skew      (N x 1) [ms]
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Updated: 23 Jan 2026
%       Updated: 260716 (add .freqPeak, the 1/f-corrected ripple frequency;
%                the Hilbert .freq is biased low by the aperiodic slope).

%% ========================================================================
%  ARGUMENTS & SETUP
%  ========================================================================

fs = ripp.info.fs;
nEvents = size(ripp.times, 1);
nSamples = length(rippSig.filt);

% Initialize Output Structure
ripp.amp = nan(nEvents, 1);
ripp.freq = nan(nEvents, 1);
ripp.freqEvent = nan(nEvents, 1);
ripp.freqPeak = nan(nEvents, 1);
ripp.peakProm = nan(nEvents, 1);
ripp.energy = nan(nEvents, 1);
ripp.dur = nan(nEvents, 1);
ripp.skew = nan(nEvents, 1);

% Whitened-peak window (fixed, symmetric about the peak)
nPeakWin = round(0.064 * fs);

% Convert Times to Samples (1-based indexing)
% We use max/min to ensure we don't index outside the signal bounds
startSamps = max(1, round(ripp.times(:, 1) * fs) + 1);
endSamps   = min(nSamples, round(ripp.times(:, 2) * fs) + 1);
peakSamps  = round(ripp.peakTime * fs) + 1;

% Ensure peaks are within bounds (sanity check)
peakSamps = max(1, min(nSamples, peakSamps));

% Fixed window (avoids duration bias)
winFxd = 0.020;
nWin = round(winFxd * fs);

%% ========================================================================
%  CALCULATE PARAMETERS
%  ========================================================================

for iEvent = 1:nEvents

    % Fixed window indices
    idxFxd = (peakSamps(iEvent) - nWin) : (peakSamps(iEvent) + nWin);
    idxFxd = idxFxd(idxFxd >= 1 & idxFxd <= nSamples);

    % Detected ripple indices
    idxDtct = startSamps(iEvent) : endSamps(iEvent);

    % Duration (ms)
    % Calculated directly from timestamps for precision
    ripp.dur(iEvent) = (ripp.times(iEvent, 2) - ripp.times(iEvent, 1)) * 1000;

    % Peak Amplitude (uV)
    % Instantaneous amplitude of the envelope at the exact peak index
    ripp.amp(iEvent) = rippSig.amp(peakSamps(iEvent));

    % Mean Frequency (Hz) - Hilbert instantaneous (biased low by 1/f)
    ripp.freq(iEvent) = mean(rippSig.freq(idxFxd), 'omitnan');
    ripp.freqEvent(iEvent) = mean(rippSig.freq(idxDtct), 'omitnan');

    % Peak Frequency (Hz) - whitened event-PSD peak (1/f-corrected)
    idxPk = (peakSamps(iEvent) - nPeakWin) : (peakSamps(iEvent) + nPeakWin);
    idxPk = idxPk(idxPk >= 1 & idxPk <= nSamples);
    [ripp.freqPeak(iEvent), ripp.peakProm(iEvent)] = ...
        ripp_freqPeak(rippSig.lfp(idxPk), fs);

    % Total Energy (uV^2)
    ripp.energy(iEvent) = sum(rippSig.filt(idxDtct) .^ 2, 'omitnan');

    % Skewness (ms)
    % Center of Mass of the Hilbert Envelope relative to the Peak Time.
    evtAmp = rippSig.amp(idxFxd);
    tRel = (idxFxd(:) - peakSamps(iEvent)) / fs; % Relative time in seconds
    ampSum = sum(evtAmp, 'omitnan');
    ripp.skew(iEvent) = (sum(tRel .* evtAmp, 'omitnan') / ampSum) * 1000;

end

end     % EOF


% =========================================================================
%  LOCAL: whitened event-PSD peak frequency
% =========================================================================
function [f0, prom] = ripp_freqPeak(seg, fs)
% Peak frequency AND prominence of one event's spectrum after removing the 1/f
% background. Fits a power law (log-log line) to the aperiodic part - the fit
% range brackets the ripple band but excludes it (30-500 Hz minus 70-260 Hz) so
% the oscillation does not pull the slope - divides it out, and returns the
% residual peak inside 70-260 Hz (f0, the honest ripple frequency; the Hilbert
% estimate rides the 1/f down) together with its height above the background
% (prom): ~1 when the event is a broadband transient with no real oscillation,
% >>1 when a narrowband ripple stands above the 1/f floor.

seg = seg(:);
prom = NaN;
if numel(seg) < 16, f0 = NaN; return; end

nf  = 512;
fAx = (0:nf/2)' * fs / nf;
ps  = abs(fft(detrend(seg) .* hann(numel(seg)), nf)) .^ 2;
ps  = ps(1:nf/2+1);

fitBand = (fAx >= 30 & fAx <= 500) & ~(fAx >= 70 & fAx <= 260) & fAx > 0;
if nnz(fitBand) < 4, f0 = NaN; return; end

pf   = polyfit(log(fAx(fitBand)), log(ps(fitBand) + eps), 1);
whit = ps ./ exp(polyval(pf, log(max(fAx, 1))));

rippBand = fAx >= 70 & fAx <= 260;
fRipp    = fAx(rippBand);
[prom, im] = max(whit(rippBand));
f0         = fRipp(im);

end