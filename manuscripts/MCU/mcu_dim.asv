function mcu_dim(varargin)
% MCU_DIM Analysis of dimensionality vs. firing rate changes.
%
%   MCU_DIM(...) performs dimensionality analysis on population activity
%   and correlates it with firing rate changes in MCU-KO and WT mice,
%   both in vivo and in vitro (MEA).
%
%   INPUTS:
%       varargin    - (param/value) Optional parameters:
%                     'type'    : (char) 'all', 'vivo', or 'vitro' {'all'}
%                     'dimMet'  : (char) Dimensionality method {'pr'}
%                     'thrVal'  : (num) Threshold for dimensionality {0.8}
%                     'binSize' : (num) Bin size for FR calc {0.1} (s)
%                     'flgPlot' : (log) Generate plots {true}
%
%   See also: DIM_CALC, MCU_TBLVIVO, MCU_TBLMEA

%   EXAMPLE CALL:
        % params.dimMet = 'pr';
        % params.thrVal = 0.8;
        % params.binSize = 0.1;
        % params.flgPlot = true;

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addParameter(p, 'type', 'all', @ischar);
addParameter(p, 'dimMet', 'pr', @ischar);
addParameter(p, 'thrVal', 0.8, @isnumeric);
addParameter(p, 'binSize', 0.1, @isnumeric);
addParameter(p, 'flgPlot', true, @islogical);
parse(p, varargin{:});

params = p.Results;

%% ========================================================================
%  MAIN EXECUTION
%  ========================================================================

switch lower(params.type)
    case 'all'
        run_vivo(params);
        run_vitro(params);
    case 'vivo'
        run_vivo(params);
    case 'vitro'
        run_vitro(params);
end

end

%% ========================================================================
%  IN VIVO ANALYSIS
%  ========================================================================

function run_vivo(params)

fprintf('Running In Vivo Analysis...\n');

% Configure varMap to include spike times
varMap = struct();
varMap.UnitType = 'units.type';
varMap.FR = 'fr.mfr';
varMap.spktimes = 'spikes.times';

% Load 
basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];
vars = {'spikes', 'fr', 'units'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);
[tbl, ~] = mcu_tblVivo('basepaths', basepaths, 'v', v, ...
    'varMap', varMap, 'flgClean', true);

% Exclude specific subject that only has baseline
tbl(tbl.Name == "lh137", :) = []; 

% Calculate Rate Changes (Log Ratio)
stats = grpstats(tbl, {'Name', 'Day'}, 'mean', 'DataVars', 'FR');
stats = stats(:, {'Name', 'Day', 'mean_FR'});
tblWide = unstack(stats, 'mean_FR', 'Day');
tblWide.Ratio = (tblWide.BAC3 ./ tblWide.BSL) * 100;
tblWide = tblWide(:, {'Name', 'Ratio'});

% Calculate Dimensionality during Baseline (BSL)
winBsl = [300, 315] * 60;
tblDim = dim_file(tbl(tbl.Day == 'BSL', :), winBsl, params);

% Join and Plot
tblCorr = join(tblDim, tblWide, 'Keys', 'Name');

if params.flgPlot
    plot_correlation(tblCorr, 'In Vivo');
end

end

%% ========================================================================
%  IN VITRO (MEA) ANALYSIS
%  ========================================================================

function run_vitro(params)

% Load MEA Table
[tbl, ~, ~] = mcu_tblMea();

% Calculate Ratio
% tbl = tbl_transform(tbl, 'flg0', true, 'verbose', true, 'varsInc', {'fr', 'frSs'});
tbl.Ratio = (tbl.frSs ./ tbl.fr) * 100;

tblStats = grpstats(tbl, 'Name', 'mean', 'DataVars', {'fr', 'frSs'});
tblStats.Ratio = (tblStats.mean_frSs ./ tblStats.mean_fr) * 100;
tblStats = tblStats(:, {'Name', 'Ratio'});

% Calculate Dimensionality
winBsl = [0, 15] * 60;
tblDim = dim_file(tbl, winBsl, params);

% Join and Plot
tblCorr = join(tblDim, tblStats, 'Keys', 'Name');

if params.flgPlot
    plot_correlation(tblCorr, 'MEA (In Vitro)');
end

end

%% ========================================================================
%  HELPER: DIMENSIONALITY CALCULATION
%  ========================================================================

function tblDim = dim_file(tbl, win, params)

% Unique files (Names) to process
fileNames = unique(tbl.Name);
nFiles = length(fileNames);
dimVal = nan(nFiles, 1);
dimGroup = cell(nFiles, 1);

tVec = win(1) : params.binSize : win(2);
nBins = length(tVec) - 1;

for iFile = 1 : nFiles

    name = fileNames(iFile);

    % Extract units for this file
    tblFile = tbl(tbl.Name == name, :);
    spktimes = tblFile.spktimes;

    % Group (Take from first row)
    dimGroup{iFile} = char(tblFile.Group(1));

    % Calculate FR Matrix
    nUnits = length(spktimes);
    frMat = nan(nUnits, nBins);

    for iUnit = 1 : nUnits
        counts = histcounts(spktimes{iUnit}, tVec);
        frMat(iUnit, :) = counts ./ params.binSize;
    end

    % Dimensionality
    dimVal(iFile) = dim_calc(frMat, 'method', params.dimMet, ...
        'thrVal', params.thrVal);

end

% Create Table
tblDim = table(categorical(dimGroup), fileNames, dimVal, ...
    'VariableNames', {'Group', 'Name', 'Val'});

end

%% ========================================================================
%  HELPER: PLOTTING
%  ========================================================================

function plot_correlation(tblCorr, titleStr)

hFig = figure('Color', 'w', 'Name', [titleStr ' Correlation']);
ax = axes(hFig);

% Scatter
scatter(ax, tblCorr.Val, tblCorr.Ratio, 60, 'filled', 'k');

% Labels
xlabel(ax, 'Dimensionality');
ylabel(ax, 'MFR Recovery');
title(ax, [titleStr ': Dimensionality vs. Rate Changes']);
grid(ax, 'on');

% Regression Line
hRef = lsline(ax);
hRef.Color = 'r';
hRef.LineWidth = 2;

% Statistics
[R, P] = corr(tblCorr.Val, tblCorr.Ratio, 'Type', 'Spearman');
strStat = sprintf('R = %.2f, p = %.3f', R, P);
text(ax, 0.05, 0.95, strStat, 'Units', 'normalized', ...
    'FontSize', 12, 'FontWeight', 'bold');

end
