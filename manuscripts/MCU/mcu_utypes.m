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

    % waveform metrices
    % swv = spkwv_metrics('basepath', basepath, 'flgSave', true, 'flgForce', true);

    % Spike timing metrics
    st = spktimes_metrics('spktimes', v(iPath).spikes.times, 'sunits', [],...
        'bins', {[0, Inf]}, 'flgForce', true, 'flgSave', true, 'flgAll', false);

end



%% ========================================================================
%  RE-CLASSIFY
%  ========================================================================

% get all files in study
basepaths = mcu_basepaths('all');

basepaths = mcu_basepaths('lh142');

% Create table of features for classification
fetTbl = utypes_features('basepaths', {basepaths{:}}, 'flgPlot', false);

% Classify
fetSelect = {'asym', 'hpk', 'tp'};
rsPrior = 0.97;
regVal = 0.01;
uTbl = utypes_classify('basepaths', {basepaths{:}}, ...
    'flgSave', false, 'fetTbl', fetTbl,...
    'fetSelect', fetSelect, 'regVal', regVal, 'rsPrior', rsPrior,...
    'flgPlot', true, 'flgSave', false);

% Plot waveforms
plot_wv('basepaths', basepaths)

% Plot classification
plot_utypes('basepaths', basepaths, 'flgSave', true)



%% ========================================================================
%  LOAD DATA TABLE
%  ========================================================================

basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];
uTbl = mcu_unitData();

% Plot scatter gui
cfg = mcu_cfg;
cfgGui.xVar = 'TP';
cfgGui.yVar = 'BLidor';
cfgGui.szVar = 'FR';
cfgGui.grpVar = 'UnitType';
cfgGui.clr = cfg.clr.unit([3 : -1 : 1], :);
cfgGui.alpha = 0.4;
hFig = plot_tblScatter(uTbl, 'cfg', cfgGui);

% Grab FR vs Time data
[tAxis, frTbl] = mcu_frTbl(basepaths, 'uTbl', uTbl, 'flgPlot', true);

% Clean up
dataTbl = frTbl;
dataTbl(dataTbl.UnitType == 'Other', :) = [];
dataTbl.UnitType = removecats(dataTbl.UnitType, 'Other');

mcu_dashboard(dataTbl, tAxis, 'Supervisor_Report.html');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Formula
frml = 'SWV ~ Group * UnitType + (1|Mouse)';

% organize for lme
varFld = 'tp';
[lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varFld', varFld);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot classification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get all files in study
basepaths = mcu_basepaths('all');
nPaths = length(basepaths);

% Get only WT basepaths
mNames = unique(get_mname(mcu_basepaths('wt')));
clear mPaths
for iMouse = 1 : length(mNames)
    mPaths(iMouse, :) = string(mcu_basepaths(mNames{iMouse}))';
end
basepaths = mPaths(:);
basepaths = [mcu_basepaths('wt_bsl'), mcu_basepaths('wt_bsl_ripp')];
basepaths = unique(basepaths);


altClassify = 3;
swvFld = 'tp';
stFld = 'royer';
hAx = plot_utypes('basepaths', basepaths, 'flgRaw', false,...
    'plotType', 'scatter', 'swvFld', swvFld, 'stFld', stFld,...
    'unitIdx', altClassify);
hFig = gcf;
xlabel(hAx, 'Trough-to-Peak (ms)')
ylabel(hAx, 'Burstiness Index')
axis tight
set(hAx, 'YScale', 'log')
ylim([0.1, 100])
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape',...
    'square', 'axHeight', 300, 'flgPos', true)

% Extract data
hSct = findobj(gca, 'Type', 'scatter');
iUnit = 2;
x = hSct(iUnit).XData;
y = hSct(iUnit).YData;
sz = hSct(iUnit).SizeData;


% Save
fname = ['UnitTypes~Scatter_Alt', num2str(altClassify)];
lme_save('hFig', hFig, 'fname', fname, 'frmt', {'mat', 'svg'},...
    'lmeData', [], 'lmeStats', [])


hAx = plot_utypes('basepaths', basepaths, 'flgRaw', false,...
    'plotType', 'wv', 'swvFld', swvFld, 'stFld', stFld,...
    'unitIdx', altClassify);
hFig = gcf;
xlabel(hAx, 'Time (ms)')
ylabel(hAx, 'Amplitude (a.u.)')
yticks(hAx, [])
ylim([-0.75, 0.1])
xlim([-0.5, 0.5])
hAx.Legend.Location = 'southwest';
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'tall', 'axHeight', 300)
% hAx.Legend.FontSize = 16;
hAx.Legend.ItemTokenSize = [14 10];

% Save
fname = ['UnitTypes~Wv_Alt', num2str(altClassify)];
lme_save('hFig', hFig, 'fname', fname, 'frmt', {'mat', 'svg'},...
    'lmeData', [], 'lmeStats', [])
















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create copy of files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get all files in study
basepaths = mcu_basepaths('all');
nPaths = length(basepaths);

fNames = {...
    '*.acceleration.mat',...
    '*.datInfo*',...
    '*.cell_metrics*',...
    '*.fr.mat',...
    '*.frEmg.mat',...
    '*.nrs',...
    '*.ripp.*',...
    '*.session.mat*',...
    '*.spikes.*',...
    '*.spktimes.mat',...
    '*.sr.mat',...
    '*.st_*',...
    '*.swv_metrics*',...
    '*.swv_raw*',...
    '*.units.mat',...
    '*.xml',...
    '*.sleep_states*',...
    '*.sleep_labelsMan*',...
    % '*.lfp',...
    % '*.sleep_*',...
    };

for iPath = 1 : nPaths
    basepath = basepaths{iPath};
    cp_basepath('fNames', fNames, 'basepath', basepath,...
        'overwrite', false, 'verbose', true)
end
