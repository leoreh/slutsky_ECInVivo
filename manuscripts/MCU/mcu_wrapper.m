
%% ========================================================================
%  PER FILE ANALYSIS
%  ========================================================================

% get all files in study
basepaths = mcu_basepaths('all');
% basepaths = unique([mcu_basepaths('wt_bsl'), mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl')]);
nFiles = length(basepaths);

% vars
vars = {'spikes', 'units', 'frNet'};
vars = {'sleep_states', 'spikes'};

% load state vars
v = basepaths2vars('basepaths', basepaths, 'vars', vars);


for iFile = 1 : nFiles

    % file
    basepath = basepaths{iFile};
    [~, basename] = fileparts(basepath);
    cd(basepath)

    % Spktimes, maybe limit to unit type and vigilance state
    spktimes = v(iFile).spikes.times;
    % uType = v(iFile).units.type;
    % spktimes = spktimes(uType == "RS");

    % FR network stats during Baseline
    % winBsl = [300, 360] * 60;
    % frNet = fr_network(spktimes, 'flgSave', false, 'winLim', winBsl);

    % winLim = [0, Inf] * 60 * 60;
    % drft = drift_file(spktimes, 'flgSave', true, 'winLim', winLim, ...
    %     'binSize', 1200, 'winSize', [], 'flgPlot', true);
    
    % Limit to state
    boutTimes = intervals(v(iFile).ss.bouts.times{4});
    spktimes = cellfun(@(x) boutTimes.restrict(x, 'flgShift', false), spktimes, 'uni', false);

    % Burst detection
    isiStart = 0.006;
    minSpks = 2;
    % isiStart = 0.05;
    % minSpks = 3;
    isiEnd = isiStart * 2;
    minIBI = isiEnd;
    minDur = 0;
    brst = brst_detect(spktimes, ...
        'minSpks', minSpks, ...
        'isiStart', isiStart, ...
        'isiEnd', isiEnd * 2, ...
        'minDur', minDur, ...
        'minIBI', minIBI, ...
        'flgForce', true, 'flgSave', true, 'flgPlot', false);

    % Burst statistics
    stats = brst_stats(brst, spktimes, 'winCalc', [], ...
        'flgSave', true);

end



%% ========================================================================
%  INSPECT BASELINE
%  ========================================================================

presets = {'frNet', 'brst'};
basepaths = [mcu_basepaths('wt_bsl'), mcu_basepaths('mcu_bsl')];
tbl = mcu_tblVivo('basepaths', basepaths, 'presets', presets);

% Transform burst metrics
tblPlot = tbl;
tblPlot = tbl_trans(tblPlot, 'varsInc', {'pBspk'}, 'logBase', 'logit', 'verbose', true);
tblPlot = tbl_trans(tblPlot, 'varsInc', {'bRoy'}, 'logBase', 10, 'verbose', true);

% Limit units
uIdx = tblPlot.unitType == 'RS';
tblPlot = tblPlot(uIdx, :);

tblGUI_bar(tblPlot, 'xVar', 'Group', 'yVar', 'funcon');

tblGUI_scatHist(tblPlot, 'xVar', 'pBspk', 'yVar', 'funcon_fish', 'grpVar', 'Group');



% LME
varRsp = 'funcon_shf';
frml = [varRsp, ' ~ Group * pBspk + (1|Name)'];
[lmeMdl, lmeStats, lmeInfo, tblMdl] = lme_analyse(tblPlot, frml);



% PLOT INTERACTION
vars = {'pBspk', 'Group'};
hFig = plot_axSize('flgFullscreen', true, 'flgPos', true);

% Partial Residuals
hAx = nexttile;
[prRes, hFig] = lme_pr(lmeMdl, vars, 'transParams', [], ...
    'hAx', hAx);

% Partial Dependence
hAx = nexttile;
vars = {'Group'};
[pdRes, hFig] = lme_pd(lmeMdl, vars, 'transParams', [], ...
    'hAx', hAx);


%% ========================================================================
%  INSPECT RECOVERY
%  ========================================================================

% Load table
basepaths = [mcu_basepaths('wt_bsl'), mcu_basepaths('mcu_bsl'), ...
    mcu_basepaths('wt_bac3'), mcu_basepaths('mcu_bac3')];
basepaths(contains(basepaths, 'lh137')) = [];
presets = {'brst'};
tblVivo = mcu_tblVivo('basepaths', basepaths, 'flgClean', true, ...
    'presets', presets);
tblVivo.Day(tblVivo.Day == "BAC_ON") = "BAC3";
tblVivo.Day = removecats(tblVivo.Day, {'BAC_ON'});

% Filter RS
tblVivo = tblVivo(tblVivo.unitType == "RS", :);
tblVivo.unitType = [];

% logit pBspk
tblTrans = tbl_trans(tblVivo, 'varsInc', {'pBspk'}, 'logBase', 'logit');
tblVivo.pBspk_trans = tblTrans.pBspk;

% Assert no zero values (instead of a pseudocount)
tblVivo = tbl_trans(tblVivo, 'flg0', true, 'verbose', true);

% Plot
tblGUI_bar(tblVivo, 'xVar', 'Group', 'yVar', 'fr');
tblGUI_scatHist(tblVivo, 'xVar', 'pBspk', 'yVar', 'fr', 'grpVar', 'Group');


%% ========================================================================
%  WIDE TABLE
%  ========================================================================

% Unstack (wide table)
tblPlot = tblVivo;
varsTbl = tblPlot.Properties.VariableNames;
isNum = cellfun(@(x) isnumeric(tblPlot.(x)) && ~iscategorical(tblPlot.(x)), varsTbl);
varsNum = varsTbl(isNum);
tblPlot = unstack(tblPlot, varsNum, {'Day'});
tblPlot.UnitID = [];

% Organize column names
varsTbl = tblPlot.Properties.VariableNames;
varsNew = regexprep(varsTbl, '_BSL$', '');
varsNew = regexprep(varsNew, '(.*)_BAC3$', 'ss_$1');
tblPlot.Properties.VariableNames = varsNew;

% Plot
tblGUI_scatHist(tblPlot, 'xVar', 'ss_frBspk', 'yVar', 'ss_frSspk', 'grpVar', 'Group');

% LME
varRsp = 'ss_frSspk';
frml = [varRsp, ' ~ Group + (1|Name)'];
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblPlot, frml, 'dist', 'gamma');



%% ========================================================================
%  COLLAPSE UNIT TABLE BY DAY
%  ========================================================================

statType = 'mean';

if strcmp(statType, 'geomean')
    fh = @stat_geomean;
else
    fh = statType;
end

varsTbl = tblPlot.Properties.VariableNames;
tblSum = groupsummary(tblPlot, {'Group', 'Name'}, fh, ...
    vartype("numeric"));
tblSum(:, "GroupCount") = [];
varsTbl(contains(varsTbl, 'File')) = [];
tblSum.Properties.VariableNames = varsTbl;

% Add Recovery Metrics

% Relative
tblSum.dBrst_rel = log((tblSum.ss_frBspk) ./ (tblSum.frBspk));
tblSum.dSngl_rel = log((tblSum.ss_frSspk) ./ (tblSum.frSspk));
tblSum.dFr = log((tblSum.ss_fr) ./ (tblSum.fr));
tblSum.dpBspk = (tblSum.ss_pBspk_trans) - (tblSum.pBspk_trans);

% Absolute
tblSum.dBrst_abs = (tblSum.ss_frBspk - tblSum.frBspk);
tblSum.dSngl_abs = (tblSum.ss_frSspk - tblSum.frSspk);
tblSum.dFr_abs = (tblSum.ss_fr - tblSum.fr);
tblSum.dFr_prct = dFr_abs * 100;


tblGUI_scatHist(tblSum, 'xVar', 'dBrst_rel', 'yVar', 'dSngl_rel', 'grpVar', 'Group');

tblSum(:, {'Group', 'dBrst_abs', 'dSngl_abs'})



%% ========================================================================
%  PSEUDOREPLICATION
%  ========================================================================

% get all files in study
basepaths = [mcu_basepaths('wt_bsl'), mcu_basepaths('mcu_bsl')];
nFiles = length(basepaths);

% vars
vars = {'spikes', 'units', 'brst', 'frNet'};

% load state vars
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

% Initialize
vPseudo = struct();
winSize = 15 * 60;

for iFile = 1 : nFiles

    % file
    basepath = basepaths{iFile};
    [~, basename] = fileparts(basepath);
    cd(basepath)
    
    % Spktimes
    spktimes = v(iFile).spikes.times;

    % Unit Filter. DIM and MCC on all units. MFR, PBSPK, DRIFT, FUNCON only on RS
    uType = v(iFile).units.type;
    uIdx = uType == "RS";

    % Calculate Dim
    frNet = fr_network(spktimes, 'flgSave', false, 'winLim', [], ...
        'winSize', winSize);
    dim = frNet.dim;
    
    % Grab funcon from pre-calculated struct
    frNet = v(iFile).frNet;
    mcc = frNet.corr.shuffle.mcc;
    funcon = mean(frNet.corr.shuffle.funcon(uIdx, :), 'omitnan');

    % FR matrix
    nWin = size(frNet.tWin, 1);
    nUnits = length(spktimes);
    winEdges = [frNet.tWin(:, 1); frNet.tWin(end, 2)];
    frMat = nan(nUnits, nWin);
    for iUnit = 1:nUnits
        frMat(iUnit, :) = histcounts(spktimes{iUnit}, winEdges, ...
            'Normalization', 'countdensity');
    end
    frMat = frMat(uIdx, :);
    mfr = mean(frMat, 1, 'omitnan')';

    % Drift (similarity)
    drft = drift_calc(frMat, 'flgPlot', false, ...
        'thrFr', 0.005, ...
        'thrLin', [], ...
        'limUnit', []);
    dtCorr = [drft.dt_corr{1}; nan];

    % Burst statistics
    brst = v(iFile).brst;
    stats = brst_stats(brst, spktimes, 'winCalc', frNet.tWin, ...
        'flgSave', false);
    pBspk = mean(stats.pBspk(uIdx, :), 1, 'omitnan')';

    % Store
    vPseudo(iFile).mfr      = mfr;
    vPseudo(iFile).mcc      = mcc;
    vPseudo(iFile).funcon   = funcon;
    vPseudo(iFile).dim      = dim;
    vPseudo(iFile).pBspk    = pBspk;
    vPseudo(iFile).dtCorr   = dtCorr;
end

% Metadata tags
tagFiles = struct();
tagFiles.Name = get_mname(basepaths);
[~, fileNames] = fileparts(basepaths);
tagFiles.File = fileNames;

% Variable Map
varMap = struct();
varMap.mfr      = 'mfr';
varMap.mcc      = 'mcc';
varMap.dim      = 'dim';
varMap.pBspk    = 'pBspk';
varMap.dtCorr   = 'dtCorr';
varMap.funcon   = 'funcon';

% Create Table
tbl = v2tbl('v', vPseudo, 'varMap', varMap, 'tagFiles', tagFiles);

% Add Group
cfg = mcu_cfg;
tbl.Group = ones(height(tbl), 1) * 1;
tbl.Group(ismember(tbl.Name, cfg.miceMCU), :) = 2;
tbl.Group = categorical(tbl.Group, [1, 2], cfg.lbl.grp);

% Reorder columns
tbl = movevars(tbl, {'Group', 'Name', 'File'}, 'Before', 1);

% Logit pBspk
tbl = tbl_trans(tbl, 'varsInc', {'pBspk'}, 'logBase', 'logit');

% Visualize
tblGUI_scatHist(tbl, 'xVar', 'mfr', 'yVar', 'dim', 'grpVar', 'Group');

% Analysis
frml = 'dim ~ (funcon + mfr + pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tbl, frml);






%% ========================================================================
%  DIM VS RCV
%  ========================================================================

% get all files in study
basepaths = [mcu_basepaths('wt_bsl'), mcu_basepaths('mcu_bsl')];
nFiles = length(basepaths);

% vars
vars = {'spikes', 'units'};

% load state vars
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

for iFile = 1 : nFiles

    % file
    basepath = basepaths{iFile};
    [~, basename] = fileparts(basepath);
    cd(basepath)
    
    % Spktimes
    spktimes = v(iFile).spikes.times;
    
    % Unit Filter. DIM and MCC on all units. MFR, PBSPK, DRIFT, FUNCON only on RS
    % uType = v(iFile).units.type;
    % uIdx = uType == "RS";
    % spktimes = spktimes(uIdx);

    nUnits = length(spktimes);

    % Unit Filter. DIM and MCC on all units. MFR, PBSPK, DRIFT, FUNCON only on RS
    uType = v(iFile).units.type;
    uIdx = uType == "RS";

    % FR network stats during Baseline
    winBsl = [300, 315] * 60;
    chunks = winBsl;
    tStart = chunks(1, 1);
    tEnd   = chunks(1, 2);

    % Create Firing Rate Matrix
    binSize = 0.1;
    t = tStart : binSize : tEnd;
    nBins = length(t) - 1;
    frMat = nan(nUnits, nBins);

    for iUnit = 1:nUnits
        st = spktimes{iUnit};
        st = st(st >= tStart & st <= tEnd);
        counts = histcounts(st, t);
        frMat(iUnit, :) = counts ./ binSize;
    end

    % Dimensionality
    dim(iFile) = dim_calc(frMat, ...
        'method', 'pr', ...
        'thrVal', 0.8, ...
        'flgFrac', true, ...
        'flgStable', true, ...
        'flgPlot', false);   

end


tblLme.dim = dim([1 : 9, 11])';

tblGUI_scatHist(tblLme, 'xVar', 'dim', 'yVar', 'Rcv', 'grpVar', 'Group');


