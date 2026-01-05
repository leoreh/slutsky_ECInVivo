function [tAxis, tblUnit] = mcu_frTbl(basepaths, varargin)

% MCU_FRTBL Create a table of firing rates aligned to perturbation
%   Combines data from multiple mice into a single global time axis.
%   Alignment Strategy:
%     1. Files 1-5: Aligned such that Perturbation Onset = t(0).
%     2. Files 6-7: Aligned such that Start of File 6 = t(144h).
%
% INPUTS:
%   basepaths (cell) - List of recording session paths.
%   flgPlot   (log)  - Results plotting. {true}
%
% OUTPUTS:
%   tAxis     (vec)  - Global time axis in Hours (0 = Perturbation).
%   tblUnit   (tbl)  - Table with metadata and 'FRt' column (aligned matrix).
%
% See also: MCU_DETECTPERT, CATFRTIME, V2TBL

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'basepaths', @iscell);
addOptional(p, 'flgPlot', true, @islogical);
parse(p, basepaths, varargin{:});

flgPlot = p.Results.flgPlot;

% Sort basepaths naturally
basepaths = natsort(basepaths);


%% ========================================================================
%  CREATE TABLE
%  ========================================================================

tblUnit = mcu_tblVivo('basepaths', basepaths);

%% ========================================================================
%  GLOBAL TIME AXIS
%  ========================================================================

% Define Global Axis Range. We create a contiguous axis that covers the min
% to max of the entire experiment, even if the middle is empty (NaNs).

binSize = 60;           % Firing rate binsize [s]
minT = -48;             % Start check 24h before pert
maxW = 146;             % End check 8 days later
dt = binSize / 3600;    % [hr]

tAxis = minT : dt : maxW;
nBinsTotal = length(tAxis);

% Initialize mat
tblUnit.FRt = nan(height(tblUnit), nBinsTotal);


%% ========================================================================
%  PROCESS BY MOUSE
%  ========================================================================

mice = natsort(unique(get_mname(basepaths)));
nMice = length(mice);
unitCount = 1;

for iMouse = 1:nMice
    mName = string(mice(iMouse));
    mPaths = basepaths(contains(basepaths, mName));

    % --- BASELINE + BACLOFEN ---
    paths = mPaths(1:5);
    frMat = cat_fr(paths);

    % Detect Perturbation (returns index relative to start of frBac)
    [pertIdx, tMouse] = mcu_detectPert(frMat, 'flgPlot', false);

    % Remove NaN samples before the perturbation. No need to adhere to
    % absolute time of day for baseline data. Update pertIdx to the shorter
    % matrix
    idxRm = all(isnan(frMat(:, tMouse < 0)), 1);
    frMat(:, idxRm) = [];
    [nUnits, nBins] = size(frMat);
    pertIdx = pertIdx - sum(idxRm);

    % Map to Global Axis: t=0 is at pertIdx
    tStart = (1 - pertIdx) * dt;

    % Find index in tAxis corresponding to peturbation onset and fill
    % matrix
    [~, idxStart] = min(abs(tAxis - tStart));
    idxEnd = min(nBinsTotal, idxStart + nBins - 1);
    tblUnit.FRt(unitCount : unitCount + nUnits - 1, idxStart : idxEnd) = frMat;
    unitCount = unitCount + nUnits;

    % --- WASHOUT ---
    paths = mPaths(6:7);
    frMat = cat_fr(paths);
    [nUnits, nBins] = size(frMat);

    % Alignment: Fixed start 4 days after perturbation onset
    tStart = 96;
    [~, idxStart] = min(abs(tAxis - tStart));
    idxEnd = min(nBinsTotal, idxStart + nBins - 1);
    tblUnit.FRt(unitCount : unitCount + nUnits - 1, idxStart : idxEnd) = frMat;
    unitCount = unitCount + nUnits;

end


% --- DENOISE ---
tblUnit.FRt = fr_denoise(tblUnit.FRt, 'flgPlot', false, 'frameLen', 60);

%% ========================================================================
%  PLOT
%  ========================================================================

if flgPlot

    tblGUI_xy(tAxis, tblUnit, 'yVar', 'FRt');

end

end
