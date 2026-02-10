%% ========================================================================
%  NOTE: PLA
%  ========================================================================

% MEA is recorded as a piecewise linear compression of time (PLA).
%   - Baseline: 20 min every hr (3x scaling) for ~4hr
%     if t_abs < 14400 | t_rec < 4800
%           t_rec = t_abs / 3
%   - Post-Perturbation: 20 min every two hr (6x scaling) for ~48hr
%     if t_abs > 14400
%           t_rec = 4800 + (t_abs - 14400) / 6
%           t_abs = 14400 + (t_rec - 4800) * 6
%
% According to spktimes (last recorded spike), recordings are ~34804s.
%   - Perturbation should occur at t_rec = 4800s = 80min
%   - Recording duration should be 14400 + 30000 * 6 (t_abs = 54hr)
%
% To calculate X hr (t_abs) in spktimes (t_rec):
%   - X = 52; 4800 + (X * 3600 - 14400) / 6
%   - winExp = [0, 33600];

%% ========================================================================
%  NOTE: SPECIFIC FILES
%  ========================================================================
% MCU-KO2 (2020.10.19_MCUKO_Bac10uM_SORTED) has a shorter recording then
% all others (Last spktime @ 33807s). PertDetect shows this stems from a
% shorter baseline (by ~16-17 min). Hence, for all recordings window length
% is 60 bins (~3hr abs)
%
% When calculating steady-state, same length (nBins) is taken for average
% stats, even though this represents x2 absolute time. Position of windows
% (bsl, trough, ss) determined per file, according to the window in
% mea_frRcv
%
% MCU-KO are perturbed less then Control (uPert: 90%  vs 95%). However, bsl
% FR is also higher (though not significant).
% 
% MCU-KO3, uID 6. Because high FR, can distort summary plots but still, do
% no exclude because interesting dynamics: extremely fast recovery
% following gradual descent. 
% 
% Inspecting BSL vs Trough FR per unit shows they are linearly correlated
% in log space, meaning the effect of BAC is a constant percentage of BSL
% FR. Indeed, perturbation depth is very weakly correlated with BSL FR.


%% ========================================================================
%  NOTE: CLEANING UNITS
%  ========================================================================
% Manual inspection of FR traces (jan 2026):
%   - No justification for limiting bsl FR.
%   - Immediately after bac (idxPert + 5 : idxPert + 10), two units reach
%     FR ~ 50 Hz. Removed in frRcv. Other units that not perturbed are
%     technically fine (captured by uPert)


%% ========================================================================
%  PER FILE ANALYSIS
%  ========================================================================

% Files
basepaths = [mcu_basepaths('mea_bac'), mcu_basepaths('mea_mcuko')];
vars = {'mea', 'fr', 'frFit', 'brst', 'rcv'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);
nFiles = length(basepaths);

% Params
isiStart = 0.05;
isiEnd = isiStart * 2;
minIbi = isiEnd * 2;
binSize = 60;
winExp = [0, 33600];

for iFile = 1 : nFiles

    % File
    basepath = basepaths{iFile};
    [~, basename] = fileparts(basepath);
    cd(basepath);

    % % Organize raw spike times
    % files = dir('*sorted*');
    % mea = mea_orgNex('fname', files.name, 'basepath', basepath, ...
    %     'flgForce', true, 'flgSave', false, 'flgPlot', true);

    % Limit to experimental window
    spktimes = v(iFile).mea.spktimes;
    lastspike = max(cellfun(@max, spktimes, 'uni', true));
    spktimes = cellfun(@(x) x(x >= winExp(1) & x <= winExp(2)), ...
        spktimes, 'UniformOutput', false);

    % Network stats
    % rcv = v(iFile).rcv;
    % winLim = rcv.info.winBsl;
    % frNet = fr_network(spktimes, 'flgSave', true, 'winLim', winLim, ...
    %     'winSize', 15 * 60);
    % drft = drift_file(spktimes, 'flgSave', false, 'winLim', winBsl, ...
    %     'binSize', 5 * 60, 'winSize', 20 * 60, 'flgPlot', true);

    % % Firing rate
    % fr = mea_frPrep(spktimes, 'binSize', binSize, ...
    %     'flgSave', true, 'flgPlot', false);

    % % FR recovery
    % fr = v(iFile).fr;
    % rcv = mea_frRcv(fr.t, fr.fr, ...
    %     'idxTrough', fr.info.idxTrough, ...
    %     'binSize', binSize, 'flgSave', true, 'flgPlot', false);

    % % FR model fit an recovery
    % frFit = mea_frFit(fr.fr, fr.t, 'FilterLen', [], ...
    %     'idxTrough', fr.info.idxTrough, ...
    %     'flgPlot', false, 'flgSave', true);
    % frFit = v(iFile).frFit;
    % rcvMdl = mea_frRcv(fr.t, frFit.frMdl, ...
    %     'idxTrough', fr.info.idxTrough, ...
    %     'binSize', binSize, 'flgSave', false, 'flgPlot', false);
    % save(fullfile(basepath, [basename, '.frRcv_mdl.mat']), 'rcvMdl', '-v7.3');

    % Burst detection
    % brst = brst_detect(spktimes, ...
    %     'minSpks', 3, ...
    %     'isiStart', isiStart, ...
    %     'isiEnd', isiEnd, ...
    %     'minDur', 0, ...
    %     'minIBI', minIbi, ...
    %     'flgForce', true, 'flgSave', true, 'flgPlot', false);

    % % Burst statistics
    % rcv = v(iFile).rcv;
    brst = v(iFile).brst;
    % % winBsl = [0, rcv.info.winBsl(2)];
    % winCalc = [rcv.info.winBsl; rcv.info.winTrough; rcv.info.winSs];
    % stats = brst_stats(brst, spktimes, 'winCalc', winCalc, ...
    %     'flgSave', true);

    % Burst temporal dynamics
    dyn = brst_dynamics(brst, spktimes, 'binSize', 60, ...
        'binSize', binSize, 'flgSave', true, 'flgPlot', false);

    % Tranfer function spikes to Ca2+
    % ca = spk2ca(spktimes, 'winCalc', [0, Inf], ...
    %     'flgPlot', true, 'flgSave', true);

end


%% ========================================================================
%  LOAD TABLE
%  ========================================================================

presets = {'time', 'steadyState', 'frNet', 'rcv', 'spktimes'};
[tbl, xVec, basepaths, v] = mcu_tblMea('presets', presets([2]));

tblLme = tbl;

% Add logit pBspk
tblTrans = tbl_trans(tblLme, 'varsInc', {'pBspk'}, 'logBase', 'logit');
tblLme.pBspk_trans = tblTrans.pBspk;


%% ========================================================================
%  PLOTS
%  ========================================================================


tblGUI_xy(xVec, tbl);

tblGUI_scatHist(tblLme, 'xVar', 'fr', 'yVar', 'pBSpk_trans', 'grpVar', 'Group');

tblGUI_bar(tbl, 'yVar', 'pBspk', 'xVar', 'Group');

tblGUI_raster(tbl, 'grpVar', 'Name', 'grpVal', 'mcu-ko2')


%% ========================================================================
%  LME - SPIKING
%  ========================================================================


% Fit
varRsp = 'frBspk';
frml = [varRsp, ' ~ Group + (1|Name)'];
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tbl, frml);

% Fit
varRsp = 'pBspk';
frml = [varRsp, ' ~ Group * fr + (1|Name)'];
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'dist', 'logit-normal');

% Plot
hFig = tblGUI_bar(tblLme, 'yVar', varRsp, 'xVar', 'Group');

% Prism
[prismMat] = tbl2prism(tblLme, 'yVar', varRsp, 'grpVar', 'Group');
mean((prismMat), 'omitnan')
mean(log10(prismMat), 'omitnan')

% Save
fname = sprintf('MEA~%s~Group', varRsp);
lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeData', [], 'lmeStats', lmeStats, 'lmeMdl', lmeMdl)


% Compare FRbsl with FRss
tblLong = stack(tblLme, {'fr', 'frSs'}, ...
    'NewDataVariableName', 'FiringRate', ...
    'IndexVariableName', 'Epoch');

varRsp = 'FiringRate';
frml = [varRsp, ' ~ Group * Epoch + (1|Name) + (1 | Name:UnitID)'];
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLong, frml, 'dist', 'gamma');


%% ========================================================================
%  COLLAPSE PER FILE
%  ========================================================================


cfg = mcu_cfg();

% Variable names
varsTbl = tbl.Properties.VariableNames;
isNum = cellfun(@(x) isnumeric(tbl.(x)) && ~iscategorical(tbl.(x)), varsTbl);
varsNum = varsTbl(isNum);

% Select Units
tblLme = tbl;
tblLme(:, "UnitID") = [];
tblLme(:, "uRcv") = [];
tblLme(:, "uPert") = [];

% Baseline Table
varsTbl = tblLme.Properties.VariableNames;
tblLme = groupsummary(tblLme, {'Name', 'Group'}, 'mean', ...
    vartype("numeric"));
tblLme(:, "GroupCount") = [];
tblLme.Properties.VariableNames = varsTbl;

tblGUI_scatHist(tblLme, 'xVar', 'dim', 'yVar', 'Rcv', 'grpVar', 'Group');















