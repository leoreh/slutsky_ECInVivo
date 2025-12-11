function [frMat, tAxis, pertIdx] = mcu_detectPert(basepaths)

% MCU_DETECTPERT Detects perturbation onset and returns aligned time axis.
%
% Usage:
%   [pertIdx, tAxis] = mcu_detectPert(basepaths)
%
% INPUTS:
%   basepaths    (cell array) Full paths to recording folders.
%
% OUTPUTS:
%   pertIdx      (double) Index of the perturbation onset relative to the
%                concatenated firing rate matrix.
%   tAxis        (double) Time axis vector (in Minutes), zero-centered at
%                the perturbation onset.
%
% See also: CATFRTIME, MCU_FRTBL

%% ========================================================================
%  GET CONCATENATED FIRING RATES
%  ========================================================================

% Concatenate data
frMat = cat_fr(basepaths);

%% ========================================================================
%  CALCULATE POPULATION MEAN
%  ========================================================================

% Calculate mean firing rate across all units (rows)
mfr = mean(frMat, 1, 'omitnan');

% Handle NaNs (e.g. gaps in recording) by filling with 0 or interpolating
mfr = fillmissing(mfr, 'knn');

%% ========================================================================
%  DETECT PERTURBATION ONSET
%  ========================================================================

% Denoise
% Note: fr_denoise requires a time vector, we provide indices here.
mfrSmooth = fr_denoise(mfr, 1:length(mfr), ...
    'flgPlot', false, 'frameLenSec', 300);

% Limit search to first n bins (or length) to find onset
searchStart = 12 * 60;
searchWin = searchStart : 96 * 60;

% Define Window (e.g., 4 hours)
win = 4 * 60;

% Look Backward: Average of previous 'win' bins
muPre = movmean(mfrSmooth(searchWin), [win 0], 'Endpoints', 'fill');

% Look Forward: Average of next 'win' bins
muPost = movmean(mfrSmooth(searchWin), [0 win], 'Endpoints', 'fill');

% Subtract: Positive peak = Biggest Drop
diffTrace = muPre - muPost;

% Find the Peak
validIdx = (win + 1) : (length(diffTrace) - win);
[~, loc] = max(diffTrace(validIdx));
pertIdx = validIdx(1) + loc - 1;
pertIdx = searchStart + pertIdx - 1;


%% ========================================================================
%  CONSTRUCT TIME AXIS
%  ========================================================================

% Get binsize from the first file to scale the time axis
v = basepaths2vars('basepaths', basepaths(1), 'vars', {'fr'});
binsize = v(1).fr.info.binsize;
dt_hr = binsize / 60 / 60;

% Create axis centered at perturbation
nBins = size(frMat, 2);
tAxis = ((1:nBins) - pertIdx) * dt_hr;

end

% EOF
