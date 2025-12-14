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
    
    fr = calc_fr(spikes.times, 'basepath', basepath,...
        'graphics', false, 'binsize', 60, 'saveVar', true,...
        'smet', 'GK', 'winBL', [0, Inf], 'winCalc', [0, Inf], 'forceA', true);

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
basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];

% Classify
fetSelect = {'Asym', 'Hpk', 'TP'};
rsPrior = 0.97;
regVal = 0.01;
tblUnit = utypes_classify('basepaths', basepaths, ...
    'fetSelect', fetSelect, 'regVal', regVal, ...
    'rsPrior', rsPrior, 'flgPlot', false);


%% ========================================================================
%  INSPECT
%  ========================================================================

% Inspect
basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];

% Grab FR vs Time data
[tAxis, tblUnit] = mcu_frTbl(basepaths, 'flgPlot', true);

% Plot classification
utypes_gui('basepaths', basepaths, 'tAxis', tAxis, 'tblUnit', tblUnit)


hFig = tblGUI_xy(tAxis, tblUnit);

% Grab to prism
idxUnits = tblUnit.UnitType == 'FS' & tblUnit.Group == 'Control';
frMat = tblUnit.FRt(idxUnits, :)';


