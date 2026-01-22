function ripp = ripp_params(rippSig, ripp)
% RIPP_PARAMS Calculates specific ripple parameters (Amp, Freq, Energy).
%
% SUMMARY:
%   Calculates "Best Practice" parameters for each ripple event:
%     1. Peak Amplitude: Instantaneous amplitude at the peak.
%     2. Mean Frequency: Average frequency over the full event dur.
%     3. Total Energy:   Sum of squared signal over the full event dur.
%
% INPUTS:
%   rippSig     - Structure with fields: .filt (filtered LFP), .amp (envelope),
%                 and .freq (instantaneous frequency).
%   ripp.times   - (N x 2) Ripple start and end times [s].
%   ripp.peakTime   - (N x 1) Ripple peak times [s].
%   fs          - (Num) Sampling frequency [Hz].
%
% OUTPUT:
%   ripp      - Structure with fields:
%                 .amp       (N x 1) [uV]
%                 .freq      (N x 1) [Hz]
%                 .energy    (N x 1) [uV^2]
%                 .dur       (N x 1) [ms]
%
% DEPENDENCIES: None.

%% ========================================================================
%  ARGUMENTS & SETUP
%  ========================================================================

fs = ripp.info.fs;
nEvents = size(ripp.times, 1);
nSamples = length(rippSig.filt);

% Initialize Output Structure
ripp.amp = nan(nEvents, 1);
ripp.freq = nan(nEvents, 1);
ripp.energy = nan(nEvents, 1);
ripp.dur = nan(nEvents, 1);

% Convert Times to Samples (1-based indexing)
% We use max/min to ensure we don't index outside the signal bounds
startSamps = max(1, round(ripp.times(:, 1) * fs) + 1);
endSamps   = min(nSamples, round(ripp.times(:, 2) * fs) + 1);
peakSamps  = round(ripp.peakTime * fs) + 1;

% Ensure peaks are within bounds (sanity check)
peakSamps = max(1, min(nSamples, peakSamps));


%% ========================================================================
%  CALCULATE PARAMETERS
%  ========================================================================

for iEvent = 1:nEvents
    
    % Duration (ms)
    % Calculated directly from timestamps for precision
    ripp.dur(iEvent) = (ripp.times(iEvent, 2) - ripp.times(iEvent, 1)) * 1000;
    
    % Peak Amplitude (uV)
    % Instantaneous amplitude of the envelope at the exact peak index
    ripp.amp(iEvent) = rippSig.amp(peakSamps(iEvent));
    
    % Define the full event dur window
    idxDur = startSamps(iEvent) : endSamps(iEvent);
    
    % Check for valid window indices
    if isempty(idxDur), continue; end
    
    % Mean Frequency (Hz)
    % Average of the instantaneous frequency across the entire event
    if isfield(rippSig, 'freq')
        ripp.freq(iEvent) = mean(rippSig.freq(idxDur), 'omitnan');
    end
    
    % Total Energy (uV^2)
    % Sum of the squared filtered signal over the dur.
    % This represents the total integrated power of the event.
    ripp.energy(iEvent) = sum(rippSig.filt(idxDur) .^ 2, 'omitnan');
    
end

end     % EOF