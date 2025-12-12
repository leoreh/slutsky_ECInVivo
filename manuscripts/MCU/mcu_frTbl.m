function [tAxis, frTbl] = mcu_frTbl(basepaths, varargin)

% MCU_FRTBL Create a table of firing rates aligned to perturbation
%   Combines data from multiple mice into a single global time axis.
%
% INPUTS:
%   basepaths (cell) - List of recording session paths.
%   uTbl      (table)- Optional. If provided, joins with the created table.
%   flgPlot   (log)  - Resules plotting.
%
% OUTPUTS:
%   tAxis     (vec)  - Global time axis in Hours (0 = Perturbation).
%   frTbl     (tbl)  - Table with metadata and 'FRt' column (aligned matrix).
%
% See also: MCU_DETECTPERT, CATFRTIME, V2TBL

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'basepaths', @iscell);
addOptional(p, 'uTbl', table(), @istable);
addOptional(p, 'flgPlot', true, @islogical);
parse(p, basepaths, varargin{:});

uTbl = p.Results.uTbl;
flgPlot = p.Results.flgPlot;

% Sort basepaths naturally
basepaths = natsort(basepaths);

%% ========================================================================
%  CREATE TABLE
%  ========================================================================

if isempty(uTbl)
    % Create table for all basepaths at once
    v = basepaths2vars('basepaths', basepaths, 'vars', {'units'});

    varMap.UnitType = 'units.type';
    tagFiles.Name = get_mname(basepaths);
    [~, fNames] = fileparts(basepaths);
    tagFiles.File = fNames;

    frTbl = v2tbl('v', v, 'varMap', varMap, 'tagFiles', tagFiles);
else
    frTbl = uTbl;
end

%% ========================================================================
%  PROCESS BY MOUSE
%  ========================================================================

mice = unique(get_mname(basepaths));
nMice = length(mice);
mData = struct();

for iMouse = 1:nMice
    mName = string(mice(iMouse));
    myPaths = basepaths(contains(basepaths, mName));

    % Perturbation Detection & Alignment
    frMat = cat_fr(myPaths);
    [~, tMouse] = mcu_detectPert(frMat, 'flgPlot', false);

    % denoise frMat
    frMat = fr_denoise(frMat, 'flgPlot', false, 'frameLen', 60);

    mData(iMouse).frMat = frMat;
    mData(iMouse).tAxis = tMouse;
    mData(iMouse).name  = mName;
end

%% ========================================================================
%  GLOBAL ALIGNMENT
%  ========================================================================

% Determine Global Time Range
minT = inf;
maxT = -inf;
dt = [];

for iMouse = 1:nMice
    t = mData(iMouse).tAxis;
    minT = min(minT, min(t));
    maxT = max(maxT, max(t));

    if isempty(dt) && length(t) > 1
        dt = t(2) - t(1);
    end
end

% Create Global Axis
tAxis = minT : dt : maxT;
nBins = length(tAxis);

% Align Data to Global Axis
nTotalUnits = height(frTbl);
frTbl.FRt = nan(nTotalUnits, nBins);

for iMouse = 1:nMice
    mName = mData(iMouse).name;
    fr = mData(iMouse).frMat;
    t = mData(iMouse).tAxis;

    % Find Start Index in Global Grid
    [~, idxStart] = min(abs(tAxis - t(1)));

    [nUnits, nTimeMouse] = size(fr);
    alignedMat = nan(nUnits, nBins);

    idxEnd = min(nBins, idxStart + nTimeMouse - 1);
    lenFill = idxEnd - idxStart + 1;
    alignedMat(:, idxStart:idxEnd) = fr(:, 1:lenFill);

    % Assign to Table
    % Direct assignment assuming unique mouse names
    mIdx = frTbl.Name == mName;
    frTbl.FRt(mIdx, :) = alignedMat;
end

%% ========================================================================
%  PLOT
%  ========================================================================

if flgPlot
    mice = unique(frTbl.Name);
    nMice = length(mice);
    cfg = mcu_cfg();

    hFig = figure('Name', 'Firing Rate vs Time', 'NumberTitle', 'off', ...
        'Position', [100 100 1000 800]);

    tabgp = uitabgroup(hFig);

    for iMouse = 1:nMice
        mName = string(mice(iMouse));
        hTab = uitab(tabgp, 'Title', mName);

        % Get data for this mouse
        mIdx = frTbl.Name == mName;
        subTbl = frTbl(mIdx, :);

        for iUnit = 1 : 2
            hAx(iMouse, iUnit) = plot_unitFR(hTab, tAxis, subTbl, cfg, iUnit);
        end

        linkaxes(hAx(iMouse, :), 'xy');
        axis tight
    end
end

end     % EoF


%% ========================================================================
%  HELPER: PLOT FR vs TIME
%  ========================================================================

function hAx = plot_unitFR(hTab, tAxis, subTbl, cfg, iUnit)

hAx = subplot(2, 1, iUnit, 'Parent', hTab);
uType = cfg.lbl.unit{iUnit};
uIdx = subTbl.UnitType == uType;
frData = subTbl.FRt(uIdx, :);
clr = cfg.clr.unit(iUnit, :);

% Plot All Traces (Light Gray)
hold(hAx, 'on');
plot(hAx, tAxis, frData', 'Color', [0.7 0.7 0.7 0.2]);

% Plot Mean (Blue)
mfr = mean(frData, 1, 'omitnan');
plot(hAx, tAxis, mfr, 'Color', clr, 'LineWidth', 2, ...
    'DisplayName', 'Mean FR');

title(hAx, sprintf('%s Units', uType), 'Interpreter', 'none');
xlabel(hAx, 'Time (Hours)');
ylabel(hAx, 'Firing Rate (Hz)');
grid(hAx, 'on');
xline(hAx, 0, '--k', 'Perturbation');
ylim(hAx, [0, ceil(max(mfr))]);
xticks(-4 * 24 : 24 : 6 * 24)

end
