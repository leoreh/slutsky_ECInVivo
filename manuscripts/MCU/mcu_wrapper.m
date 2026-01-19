
%% ========================================================================
%  PER FILE ANALYSIS
%  ========================================================================

% get all files in study
basepaths = mcu_basepaths('all');
% basepaths = [mcu_basepaths('wt_bsl'), mcu_basepaths('mcu_bsl')];
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

    % Spktimes, maybe limit to unit type and vigilance state
    spktimes = v(iFile).spikes.times;
    % uType = v(iFile).units.type;
    % spktimes = spktimes(uType == "RS");

    % FR network stats during Baseline
    winBsl = [300, 315] * 60;
    frNet = fr_network(spktimes, 'flgSave', false, 'winLim', winBsl);

    % winLim = [0, Inf] * 60 * 60;
    % drft = drift_file(spktimes, 'flgSave', true, 'winLim', winLim, ...
    %     'binSize', 1200, 'winSize', [], 'flgPlot', true);

    % % Burst detection
    % isiVal = 0.02;
    % brst = brst_detect(spktimes, ...
    %     'minSpks', 3, ...
    %     'isiStart', isiVal, ...
    %     'isiEnd', isiVal * 2, ...
    %     'minDur', 0.005, ...
    %     'minIBI', 0.05, ...
    %     'flgForce', true, 'flgSave', true, 'flgPlot', false);

    % % Burst statistics
    % stats = brst_stats(brst, spktimes, 'winCalc', [], ...
    %     'flgSave', true);

end

%% ========================================================================
%  INSPECT BASELINE
%  ========================================================================

presets = {'frNet', 'brst'};
basepaths = [mcu_basepaths('wt_bsl'), mcu_basepaths('mcu_bsl')];
tbl = mcu_tblVivo('basepaths', basepaths, 'presets', presets);

% Logit pBspk
tbl = tbl_trans(tbl, 'varsInc', {'pBspk'}, 'logBase', 'logit');

% Limit units
tblPlot = tbl;
uIdx = tblPlot.UnitType == 'RS';
tblPlot = tblPlot(uIdx, :);

tblGUI_bar(tblPlot, 'xVar', 'Group', 'yVar', 'funcon');

tblGUI_scatHist(tblPlot, 'xVar', 'pBspk', 'yVar', 'funcon_fish', 'grpVar', 'Group');



% LME
varRsp = 'funcon_shf';
varRsp = 'funcon_fish';
frml = [varRsp, ' ~ Group  + (1|Name)'];
[lmeMdl, lmeStats, lmeInfo, tblMdl] = lme_analyse(tblPlot, frml);



%% ========================================================================
%  INSPECT RECOVERY
%  ========================================================================

% Grab only BSL and BAC3 (1st and 5th recording day)
basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];
fileMouse = findgroups(get_mname(basepaths));
[~, idxMouse] = unique(fileMouse, 'stable');
idxMouse = sort([idxMouse; idxMouse + 4]);
basepaths = basepaths(idxMouse);

presets = {'frNet', 'brst', 'prc'};
[tbl, basepaths, ~] = mcu_tblVivo('basepaths', basepaths, 'presets', presets);
tbl.Day(tbl.Day == "BAC_ON") = "BAC3";

tblPlot = tbl;
tblPlot = tbl_trans(tblPlot, 'flg0', true);

tblGUI_bar(tblPlot, 'xVar', 'Group', 'yVar', 'conn');
tblGUI_scatHist(tblPlot, 'xVar', 'fr', 'yVar', 'pBspk', 'grpVar', 'Group');



%% ========================================================================
%  COLLAPSE UNIT TABLE BY DAY
%  ========================================================================

% Variable names
varsTbl = tbl.Properties.VariableNames;
isNum = cellfun(@(x) isnumeric(tbl.(x)) && ~iscategorical(tbl.(x)), varsTbl);
varsNum = varsTbl(isNum);

% Select Units
tblRs = tbl(tbl.UnitType == "RS", :);
tblRs(:, "UnitID") = [];
tblRs(:, "File") = [];
tblRs(:, "UnitType") = [];

% Baseline Table
tblBsl = tblRs(tblRs.Day == "BSL", :);
varsTbl = tblBsl.Properties.VariableNames;
tblBsl = groupsummary(tblBsl, {'Group', 'Name', 'Day'}, 'mean', ...
    vartype("numeric"));
tblBsl(:, "GroupCount") = [];
tblBsl.Properties.VariableNames = varsTbl;

% SteadyState Table
tblSs = tblRs(tblRs.Day == "BAC3", :);
tblSs = groupsummary(tblSs, {'Group', 'Name', 'Day'}, 'mean', ...
    vartype("numeric"));
tblSs(:, "GroupCount") = [];
tblSs.Properties.VariableNames = varsTbl;

% FR Recovery
tblLme = tblBsl;
tblLme.Rcv = (tblSs.FR ./ tblBsl.FR) * 100;


tblGUI_scatHist(tblLme, 'xVar', 'dim', 'yVar', 'Rcv', 'grpVar', 'Group');



%% ========================================================================
%  PSEUDOREPLICATION
%  ========================================================================

% get all files in study
basepaths = [mcu_basepaths('wt_bsl'), mcu_basepaths('mcu_bsl')];
nFiles = length(basepaths);

% vars
vars = {'spikes', 'units', 'brst'};

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

    % FR network stats
    frNet = fr_network(spktimes, 'flgSave', true, 'winLim', [], ...
        'winSize', winSize);
    mcc = frNet.corr.shuffle.mcc;
    funcon = mean(frNet.corr.shuffle.funcon(uIdx, :), 'omitnan');
    dim = frNet.dim;

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
    iFile
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

