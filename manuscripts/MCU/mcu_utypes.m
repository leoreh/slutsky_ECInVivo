% mcu_cellCalss


%% ========================================================================
%  RE-ANALYZE
%  ========================================================================

% get all files in study
basepaths = mcu_basepaths('all');
nPaths = length(basepaths);

% vars
vars = {'spikes'};

% load state vars
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

for iPath = 1 : nPaths

    % files
    basepath = basepaths{iPath};
    [~, basename] = fileparts(basepath);
    cd(basepath)
    
    fr = spk_rate(spikes.times, 'basepath', basepath,...
        'binsize', 60, 'flgSave', true,...
        'smet', 'GK', 'winBL', [0, Inf], 'winCalc', [0, Inf]);

    % waveform metrices
    % swv = spkwv_metrics('basepath', basepath, 'flgSave', true, 'flgForce', true);

    % Spike timing metrics
    % st = spktimes_metrics('spktimes', v(iPath).spikes.times, 'sunits', [],...
    %     'bins', {[0, Inf]}, 'flgForce', true, 'flgSave', true, 'flgAll', false);

    


end



%% ========================================================================
%  RE-CLASSIFY
%  ========================================================================

% get all files in study
basepaths = mcu_basepaths('all');
basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];

% Classify
fetSelect = {'Asym', 'Hpk', 'TP'};
rsPrior = 0.97;
regVal = 0.01;
tblUnit = utypes_classify('basepaths', basepaths(3), ...
    'fetSelect', fetSelect, 'regVal', regVal, ...
    'rsPrior', rsPrior, 'flgPlot', true);


%% ========================================================================
%  INSPECT
%  ========================================================================

basepaths = [mcu_basepaths('wt_bsl_ripp')];
[tAxis, tblUnit] = mcu_frTbl(basepaths(3), 'flgPlot', false);

fileID = 9;
tblUnit = mcu_tblVivo('basepaths', basepaths(fileID), 'presets', {'swv', 'rippSpks', 'burst'});
utypes_gui('basepaths', basepaths(fileID), 'tblUnit', tblUnit)


% Inspect
basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];

% Grab FR vs Time data
[tAxis, tblUnit] = mcu_frTbl(basepaths, 'flgPlot', false);

% Plot classification
utypes_gui('basepaths', basepaths, 'tAxis', tAxis, 'tblUnit', tblUnit)


hFig = guiTbl_xy(tAxis, tblUnit);

% Grab to prism
idxUnits = tblUnit.UnitType == 'FS' & tblUnit.genotype == 'Control';
frMat = tblUnit.FRt(idxUnits, :)';


%% ========================================================================
%  AUTOCORRELOGRAM (ACG) VISUALIZATION
%  ========================================================================

% Load unit table with narrow ACG traces.
%   acgNarrow (nunits x 201): auto-correlogram at 0.5 ms resolution,
%   computed over the full recording duration (bins = [0 Inf]).
%   xAcg (1 x 201): lag axis in milliseconds, centered at 0.
basepaths = [mcu_basepaths('wt_bsl'), mcu_basepaths('mcu_bsl')];

[tblAcg, ~, ~, xAcg] = mcu_tblVivo('basepaths', basepaths, 'presets', {'acg'}, 'flgClean', true);

% Interactive viewer: tiles = unit type, colors = genotype.
% Uses 'Spread' dispersion with arithmetic mean + SEM by default.
hFig = guiTbl_xy(xAcg.narrow, tblAcg, ...
    'yVar',    'acgNarrow', ...
    'tileVar', 'unitType', ...
    'grpVar',  'genotype', ...
    'xLbl',    'Lag [ms]');


tblAcg.acgNarrow(tblAcg.genotype == 'Control')

prismTbl = groupsummary(tblAcg, "Group", {'mean', 'std'}, 'acgNarrow');

% 1. Extract Data for Control (Row 1)
mean_ctrl = prismTbl.mean_acgNarrow(1, :)'; % Transpose to column
std_ctrl  = prismTbl.std_acgNarrow(1, :)';
n_ctrl    = ones(size(mean_ctrl)) * prismTbl.GroupCount(1); % Fill column with N

% 2. Extract Data for MCU-KO (Row 2)
mean_ko = prismTbl.mean_acgNarrow(2, :)';
std_ko  = prismTbl.std_acgNarrow(2, :)';
n_ko    = ones(size(mean_ko)) * prismTbl.GroupCount(2);

% 3. Combine into the 6-column Prism format
prism_data = [mean_ctrl, std_ctrl, n_ctrl, mean_ko, std_ko, n_ko];