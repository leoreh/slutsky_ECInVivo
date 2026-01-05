

clear params
params.winLim = [0 Inf];        % Analysis window [s]
params.binSize = 0.001;         % 1ms bins
params.winStpr = 1.0;           % 1s window
params.nShuffles = 100;         % Number of shuffles
params.spkLim = Inf;

basepaths = [mcu_basepaths('wt_bsl'), mcu_basepaths('mcu_bsl')];
nFiles = length(basepaths);
vars = {'spikes', 'st_metrics', 'fr', 'units', 'spktimes', 'session', 'sleep_states'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

idxState = 1;
recDur = 60 * 60;

for iFile = 1 : nFiles

    % File
    basepath = basepaths{iFile};
    cd(basepath);

    % State
    bouts = v(iFile).ss.bouts.times{idxState};
    idxBouts = bouts(:, 1) > 6 * 60 * 60;
    bouts = bouts(idxBouts, :);

    fs = v(iFile).session.extracellular.sr;
    mutimes = cellfun(@(x) x / fs, v(iFile).spktimes, 'uni', false);
    spktimes = v(iFile).spikes.times';

    % Restrict to bouts
    spktimes = spktimes_clip(spktimes, bouts, recDur);
    mutimes = spktimes_clip(mutimes, bouts, recDur);

    tic
    prc = prc_calc(v(iFile).spikes.times, ...
        params, ...
        'mutimes', mutimes, ...
        'flgSave', true);
    toc

    prc_plot(prc, 'basepath', basepath, 'flgSave', true);

end


%  ========================================================================
% LME

% Load table
basepaths = [mcu_basepaths('wt_bsl'), mcu_basepaths('mcu_bsl')];
cfg = mcu_cfg();
tblUnit = mcu_tblVivo('basepaths', basepaths);

% Remove bad units
tblUnit(tblUnit.UnitType == 'Other', :) = [];
tblUnit.UnitType = removecats(tblUnit.UnitType, {'Other'});
tblUnit(tblUnit.UnitType == 'FS', :) = [];
tblUnit.UnitType = removecats(tblUnit.UnitType, {'FS'});
tblLme = tblUnit;

% Select Params
varRsp = 'PRC';
frml = [varRsp, ' ~ Group + (1|Name)'];
clear cfgLme
cfgLme.dist = 'Normal';

% Fit
[lmeStats, lmeMdl] = lme_analyse(tblLme, frml, cfgLme);

sum(isnan(tblLme.PRC))

% Plot
hFig = tblGUI_bar(tblLme, 'yVar', varRsp, 'xVar', 'Group');

% Prism
[prismMat] = tbl2prism(tblLme, 'yVar', varRsp, 'grpVar', 'Group');







%% ========================================================================
%  MEA
%  ========================================================================

clear params
recDur = 70 * 60;
params.winLim = [0 recDur];         % Analysis window [s]
params.binSize = 0.001;             % 1ms bins
params.winStpr = 1.0;               % 1s window
params.nShuffles = 50;              % Number of shuffles
params.spkLim = Inf;

basepaths = [mcu_basepaths('mea_bac'), mcu_basepaths('mea_mcuko')];

nFiles = length(basepaths);
vars = {'mea'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);


for iFile = 5 : nFiles

    % File
    basepath = basepaths{iFile};
    cd(basepath);

    tic
    prc = prc_calc(v(iFile).mea.spktimes, ...
        params, ...
        'flgSave', true);
    toc

    prc_plot(prc, 'basepath', basepath, 'flgSave', true);

end




[tbl, xVec, ~, v] = mcu_tblMea();

% Select Params
varRsp = 'prc';
frml = [varRsp, ' ~ Group + (1|Name)'];

% Configuration
clear cfgLme
cfgLme.dist = 'Normal';

% Fit
[lmeStats, lmeMdl] = lme_analyse(tbl, frml, cfgLme);

% Plot
hFig = tblGUI_bar(tbl, 'yVar', varRsp, 'xVar', 'Group');

% Prism
prismMat = tbl2prism(tbl, 'yVar', 'prc', 'grpVar', 'Group');
