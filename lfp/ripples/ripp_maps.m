function rippMaps = ripp_maps(rippSig, peakTime, fs, varargin)
% RIPP_MAPS Generates Peri-Event Time Histograms (PETH) for LFP signals.
%
%   rippMaps = RIPP_MAPS(rippSig, peakTime, fs, varargin)
%
%   SUMMARY:
%       Extracts windowed signal traces around ripple peaks.
%       Uses high-performance vectorized indexing to extract all events simultaneously.
%       Handles edge cases (start/end of recording) by padding with NaNs.
%
%   INPUTS:
%       rippSig     - (Struct) Structure with signal fields (.lfp, .filt, .amp, etc.).
%                              All vector fields matching LFP length will be mapped.
%       peakTime    - (Vec)    [N x 1] Times of ripple peaks [s].
%       fs          - (Num)    Sampling frequency [Hz].
%       varargin    - Parameter/Value pairs:
%           'mapDur'   - (Vec)  Window size [pre post] in seconds. (Default: [-0.05 0.05]).
%           'basepath' - (Char) Save location. (Default: pwd).
%           'flgSave'  - (Log)  Save output to .mat file? (Default: true).
%
%   OUTPUTS:
%       rippMaps    - (Struct) Contains [N_events x N_samples] matrix for each field.
%                       .tstamps - Relative time vector for the window.
%                       .lfp, .filt, .amp, etc.
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Updated: 23 Jan 2026

% =========================================================================
%  ARGUMENTS
% =========================================================================
p = inputParser;
addRequired(p, 'rippSig', @isstruct);
addRequired(p, 'peakTime', @isnumeric);
addRequired(p, 'fs', @isnumeric);
addParameter(p, 'mapDur', [-0.05 0.05], @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', true, @islogical);
parse(p, rippSig, peakTime, fs, varargin{:});

mapDur = p.Results.mapDur;
basepath = p.Results.basepath;
flgSave = p.Results.flgSave;

% =========================================================================
%  VECTORIZED EXTRACTION
% =========================================================================
[~, basename] = fileparts(basepath);
mapFile = fullfile(basepath, [basename, '.rippMaps.mat']);

% Create a relative window vector
winSamps = round(mapDur(1)*fs) : round(mapDur(2)*fs);

rippMaps = struct();
rippMaps.tstamps = winSamps / fs;
nSamps = length(rippSig.lfp);

% Convert peak times to samples
peakSamps = round(peakTime * fs) + 1;

% 4. Broadcast to create Full Index Matrix
% Creates indices for EVERY sample of EVERY ripple in one step: [nEvents x nWindowSize]
idxMat = peakSamps + winSamps;

% Handle Edge Cases (Ripples at very start/end of recording)
% We set invalid indices to 1 temporarily and then NaN them out later
validMask = (idxMat >= 1) & (idxMat <= nSamps);
idxMat(~validMask) = 1;

% 5. Extract Maps per Field
fields = fieldnames(rippSig);
for iFld = 1:length(fields)
    fn = fields{iFld};
    signal = rippSig.(fn);

    % Skip if signal dimension doesn't match LFP (e.g. metadata fields)
    if ~isvector(signal) || length(signal) ~= nSamps
        continue;
    end

    % Extract
    rawMap = signal(idxMat);

    % Apply Edge Correction (NaN out invalid indices)
    if ~all(validMask(:))
        % Cast to double to support NaNs if it was int
        rawMap = double(rawMap);
        rawMap(~validMask) = NaN;
    end

    rippMaps.(fn) = rawMap;
end

% =========================================================================
%  SAVE
% =========================================================================
if flgSave
    save(mapFile, 'rippMaps', '-v7.3');
end

end