
%  MCU_SPIKES - Analyze in vivo MCU spike data

%% ========================================================================
%  LOAD_DATA 
%  ========================================================================

basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];
cfg = mcu_cfg();

% Load table
tblUnit = mcu_unitTbl('basepaths', basepaths);

% Remove bad units
tblUnit(tblUnit.UnitType == 'Other', :) = [];
tblUnit.UnitType = removecats(tblUnit.UnitType, 'Other');

% Remove bac on / off
tblUnit(tblUnit.Day == 'BAC_ON', :) = [];
tblUnit(tblUnit.Day == 'BAC_OFF', :) = [];
tblUnit.Day = removecats(tblUnit.Day, {'BAC_ON', 'BAC_OFF'});

% Assert no zero values
tblLme = tbl_transform(tblUnit, 'flg0', true, 'verbose', true);


%% ========================================================================
%  BACLOFEN (y ~ Group * Day + (Day|Name))
%  ========================================================================

% Select Params
unitType = 'RS';
varRsp = 'FR';
varRsp = 'BRoy';

% Configuration
clear cfgLme
cfgLme.dist = 'Gamma';
cfgLme.contrasts = 1 : 19;

% Select data
tblLme = tblLme(tblLme.UnitType == unitType, :);

% run lme
frml = [varRsp, ' ~ Group * Day + (Day|Name)'];
[lmeStats, lmeMdl] = lme_analyse(tblLme, frml, cfgLme);

% Plot
hFig = tblGUI_bar(tblLme, 'yVar', varRsp, 'xVar', 'Day', 'grpVar', 'Group');


% Check best model
statsPark = lme_parkTest(tblLme, frml);
statsDist = lme_compareDists(tblLme, frml);


% Prism
iGrp = 1;
[prismMat] = tbl2prism(tblLme(tblLme.Group == cfg.lbl.grp{iGrp}, :), ...
    'yVar', varRsp, 'grpVar', 'Day');



%% ========================================================================
%  BASELINE (y ~ Group * UnitType + (1|Name))
%  ========================================================================


% Select Params
varRsp = 'FR';
varRsp = 'BRoy';

% Select data
tblLme = tblUnit(tblUnit.Day == 'BSL', :);
tblLme = tblLme(tblLme.UnitType == 'RS', :);

% Check best model
statsPark = lme_parkTest(tblLme, frml);
statsDist = lme_compareDists(tblLme, frml);

% Configuration
clear cfgLme
cfgLme.dist = 'Gamma';

% Fit
frml = [varRsp, ' ~ Group + (1|Name)'];
[lmeStats, lmeMdl] = lme_analyse(tblLme, frml, cfgLme);

% Plot
hFig = tblGUI_bar(tblLme, 'yVar', varRsp, 'xVar', 'Group');

% Prism
[prismMat] = tbl2prism(tblLme, 'yVar', varRsp);





%% ========================================================================
%  FR VS TIME
%  ========================================================================

% Files
basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];

% Grab FR vs Time data
[tAxis, tblUnit] = mcu_frTbl(basepaths, 'flgPlot', false);

% Plot
hFig = tblGUI_xy(tAxis, tblUnit);

% Grab to prism
idxUnits = tblUnit.UnitType == 'FS' & tblUnit.Group == 'Control';
frMat = tblUnit.FRt(idxUnits, :)';

prismMat = [mean(frMat, 2, 'omitnan'), ...
    std(frMat, [], 2, 'omitnan'), ...
    sum(~isnan(frMat), 2, 'omitnan')];



%% ========================================================================
%  REPRESENTATIVE RASTER
%  ========================================================================

% Files
basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];
basepaths = natsort(basepaths);
idxFiles = [1, 5, 29, 33];

% Config
cfg = mcu_cfg();
winPlot = [0, 60];
LineFormat.Color = [0 0 0];
LineFormat.LineWidth = 0.35;

% Load
vars = {'spikes', 'units'};
v = basepaths2vars('basepaths', basepaths(idxFiles), 'vars', vars);

% Plot
% close all
for iFile = 1 : length(idxFiles)

% Transpose
spktimes = cellfun(@(x) x', ...
    v(iFile).spikes.times, 'uni', false)';

% Select only RS
spktimes = spktimes(v(iFile).units.type == 'RS');

% Plot
[hFig, hAx] = plot_axSize('szOnly', false,...
    'axWidth', 600, 'axHeight', 300, 'flgPos', true);

[~, hPlt] = plot_raster(spktimes, 'PlotType', 'vertline', ...
    'lineHeight', 0.7, ...
    'lineWidth', 0.35, ...
    'hAx', hAx, ...
    'clr', [0 0 0], ...
    'xLim', winPlot, ...
    'spkDur', 0.0005);

xlim([10 15])
xlabel('Time (s)')
ylabel('Unit No.')
set(gca, 'YDir', 'normal');

[hFig, hAx] = plot_axSize('hFig', hFig, 'szOnly', false,...
    'axWidth', 600, 'axHeight', 300, 'flgPos', true);

end