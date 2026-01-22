function maps = ripp_maps(rippSig, peakTime, fs, varargin)
% RIPP_MAPS Generates Peri-Event Time Histograms (maps) for continuous signals.
%
% SUMMARY:
%   Extracts windowed signal traces around ripple peaks using fast 
%   vectorized indexing.
%
% INPUT:
%   rippSig     - Struct with fields (.lfp, .filt, .amp, etc.)
%   peakTime    - [N x 1] Times of ripple peaks [s]
%   fs          - Sampling rate [Hz]
%   varargin    - 'win' (default [-0.1 0.1]), 'basepath', 'flgSave'
%
% OUTPUT:
%   maps        - Struct with fields corresponding to rippSig fields.
%                 Each field is [nEvents x nSamples].

% =========================================================================
%  ARGUMENTS
% =========================================================================
p = inputParser;
addRequired(p, 'rippSig', @isstruct);
addRequired(p, 'peakTime', @isnumeric);
addRequired(p, 'fs', @isnumeric);
addParameter(p, 'win', [-0.1 0.1], @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', true, @islogical);
parse(p, rippSig, peakTime, fs, varargin{:});

win = p.Results.win;
basepath = p.Results.basepath;
flgSave = p.Results.flgSave;

% =========================================================================
%  VECTORIZED EXTRACTION
% =========================================================================
[~, basename] = fileparts(basepath);
mapFile = fullfile(basepath, [basename, '.rippMaps.mat']);

maps = struct();
maps.tstamps = linspace(win(1), win(2), round(diff(win)*fs)+1);
nSamps = length(rippSig.lfp);
nEvents = length(peakTime);

% Convert peak times to samples
peakSamps = round(peakTime * fs) + 1;

% Create a relative window vector
winSamps = round(win(1)*fs) : round(win(2)*fs);

% Broadcast to create a full matrix of indices [nEvents x nWindowSize]
% This creates the indices for EVERY sample of EVERY ripple in one step.
idxMat = peakSamps + winSamps; 

% Handle Edge Cases (Ripples at very start/end of recording)
% We set invalid indices to 1 temporarily and then NaN them out later
validMask = (idxMat >= 1) & (idxMat <= nSamps);
idxMat(~validMask) = 1; 

% Extract Maps
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
    
    maps.(fn) = rawMap;
end

% =========================================================================
%  SAVE
% =========================================================================
if flgSave
    save(mapFile, 'maps', '-v7.3');
    fprintf('Saved heavy maps to: %s\n', mapFile);
end

end