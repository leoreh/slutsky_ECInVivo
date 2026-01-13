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
%     FR ~ 50 Hz. Removed in frRcv.Other units that not perturbed are
%     technically fine (captured by uPert)


%% ========================================================================
%  PER FILE ANALYSIS
%  ========================================================================

% Files
basepaths = [mcu_basepaths('mea_bac'), mcu_basepaths('mea_mcuko')];
vars = {'mea', 'fr', 'frFit'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);
nFiles = length(basepaths);

% Params
isiVal = 0.05;
binSize = 60;
winExp = [0, 33600];

close all
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

    % % Network stats
    % frNet = fr_network(spktimes, 'flgSave', true, 'winLim', [0, 15] * 60);
    % drft = drift_file(spktimes, 'flgSave', false, 'winLim', winBsl, ...
    %     'binSize', 5 * 60, 'winSize', 20 * 60, 'flgPlot', true);

    % % Firing rate
    % fr = mea_frPrep(spktimes, 'binSize', binSize, ...
    %     'flgSave', true, 'flgPlot', false);

    % FR recovery
    fr = v(iFile).fr;
    rcv = mea_frRcv(fr.t, fr.fr, ...
        'idxTrough', fr.info.idxTrough, ...
        'binSize', binSize, 'flgSave', true, 'flgPlot', true);

    % FR model fit an recovery
    frFit = mea_frFit(fr.fr, fr.t, 'FilterLen', [], ...
        'idxTrough', fr.info.idxTrough, ...
        'flgPlot', false, 'flgSave', true);
    % frFit = v(iFile).frFit;
    rcvMdl = mea_frRcv(fr.t, frFit.frMdl, ...
        'idxTrough', fr.info.idxTrough, ...
        'binSize', binSize, 'flgSave', false, 'flgPlot', false);
    save(fullfile(basepath, [basename, '.frRcv_mdl.mat']), 'rcvMdl', '-v7.3');

    % Burst detection
    brst = brst_detect(spktimes, ...
        'minSpks', 3, ...
        'isiStart', isiVal, ...
        'isiEnd', isiVal * 2, ...
        'minDur', 0.005, ...
        'minIBI', 0.1, ...
        'flgForce', true, 'flgSave', true, 'flgPlot', false);

    % Burst temporal dynamics
    dyn = brst_dynamics(brst, spktimes, 'binSize', 60, 'ksd', 300, ...
        'binSize', binSize, 'flgSave', true, 'flgPlot', false);

    % Burst statistics
    winCalc = [rcv.info.winBsl; rcv.info.winTrough; rcv.info.winSs];
    stats = brst_stats(brst, spktimes, 'winCalc', winCalc, ...
        'flgSave', true);

    % Tranfer function spikes to Ca2+
    % ca = spk2ca(spktimes, 'winCalc', [0, Inf], ...
    %     'flgPlot', true, 'flgSave', true);

end


%% ========================================================================
%  LOAD TABLE
%  ========================================================================

presets = {'time', 'steadyState', 'frNet', 'rcv', 'spktimes'};
[tbl, xVec, basepaths, v] = mcu_tblMea('presets', presets([1, 3 : 4]));


%% ========================================================================
%  PLOTS
%  ========================================================================

frMin = 1 / (20 * 60);
tblLme = tbl;
tblLme.censMask = tbl.frTrough < frMin + 1e-9;
frml = 'frTrough ~ fr';

[frTroughRcv, ~] = mea_lmCens(tblLme, frml, 'censVar', 'censMask');



tblGUI_xy(xVec, tbl);

tblGUI_scatHist(tbl, 'xVar', 'pBspk', 'yVar', 'pertDepth', 'grpVar', 'Group');

tblGUI_bar(tbl, 'yVar', 'pBspk', 'xVar', 'Group');

tblGUI_raster(tbl, 'grpVar', 'Name', 'grpVal', 'mcu-ko2')



sort(tbl.frTrough)
find(tbl.frTrough > 10)

uIdx = 6;
figure
plot(xVec, tbl.t_fr(uIdx, :))

%% ========================================================================
%  LME - SPIKING
%  ========================================================================

tblLme = tbl;

% Fit
varRsp = 'frTrough';
frml = [varRsp, ' ~ Group + (1|Name)'];
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tbl, frml);

% Plot
hFig = tblGUI_bar(tblLme, 'yVar', varRsp, 'xVar', 'Group');

% Prism
[prismMat] = tbl2prism(tblLme, 'yVar', varRsp, 'grpVar', 'Group');
tblLme.fr(tblLme.Group == 'MCU-KO')

% Save
fname = sprintf('MEA~%s~Group', varRsp);
lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeData', [], 'lmeStats', lmeStats, 'lmeMdl', lmeMdl)





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

% Baseline Table
varsTbl = tblLme.Properties.VariableNames;
tblLme = groupsummary(tblLme, {'Name', 'Group'}, 'mean', ...
    vartype("numeric"));
tblLme(:, "GroupCount") = [];
tblLme.Properties.VariableNames = varsTbl;

tblGUI_scatHist(tblLme, 'xVar', 'dim', 'yVar', 'Rcv', 'grpVar', 'Group');





%% ========================================================================
%  INSPECT BURST PARAMS
%  ========================================================================

basepaths = [mcu_basepaths('mea_bac')];
vars = {'mea'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

mea = catfields([v(:).mea], 2);
spktimes = mea.spktimes;

% ISI VALLEY AS THRESHOLD
isiVal = brst_isiValley(spktimes, 'nSpks', 3);




%% ========================================================================
%  PLOT CORRELATIONS
%  ========================================================================


% -------------------------------------------------------------------------
% SELECTIVE (updated sep-25)
lData = tbl_trans(lmeData, 'varsInc', {'frBsl', 'BSpks'}, 'flgZ', false,...
    'skewThr', 0.1, 'varsGrp', {'Group'}, 'logBase', 10);
clr = mcu_cfg();
clr = clr.clr;

% Remove units that didn't recover from rcvTime
idxUnits = lmeData.uRcv;
lData.rcvTime(~idxUnits) = nan;

% Define variables now as rows and columns
varsRow = {'spkDfct', 'rcvWork', 'rcvTime'};
varsCol = {'BSpks', 'frBsl', 'pertDepth'};

% Pretty labels in the order [cols, rows] to rename properties
varsInc = [varsCol, varsRow];
clear txtLbls
txtLbls{1} = 'P(S\inB)';
txtLbls{2} = 'Firing Rate (Hz)';
txtLbls{3} = 'Perturbation Depth';
txtLbls{4} = 'Spike Deficit';
txtLbls{5} = 'Recovery Ratio';
txtLbls{6} = 'Recovery Time (Hr)';

% Restrict to included variables and rename
lData = lData(:, ['Group', varsInc]);
for iVar = 1 : length(varsInc)
    lData.Properties.VariableNames(2 : end) = txtLbls;
end

% Plot with specified rows and columns
[hFig, ~, hGrid] = plot_corrHist(lData, 'varsRow', txtLbls(4:6), ...
    'varsCol', txtLbls(1:3), 'grpIdx', 'Group', 'clrGrp', clr.grp, 'thrOut', 100);
plot_axSize('hFig', hFig, 'szOnly', false,...
    'axWidth', 1200, 'axHeight', 600, 'flgPos', true);

% Update BSpks x-tick labels on the rcvTime row (bottom-left scatter)
hAx = hGrid{4, 1};
xticks(hAx, [-2, -1, 0])
xticklabels(hAx, {'0.01', '0.1', '1'});

hAx = hGrid{4, 2};
xticks(hAx, [-2, -1, 0, 0.477])
xticklabels(hAx, {'0.01', '0.1', '1', '3'});

[nRow, nCol] = size(hGrid);
clear xLim yLim
xLim = cell(nCol, 1);
xLim{1} = [-2, 0];      % BSpks
xLim{2} = [-1, 0.5];    % frBsl (log scale representation)
xLim{3} = [5, 12];      % pertDepth
yLim = cell(nRow - 1, 1);
yLim{1} = [-4, 5];      % spkDfct
yLim{2} = [0, 1.5];   % rcvWork
yLim{3} = [-15, 40];      % rcvTime (hr)

for iRow = 1 : nRow
    for iCol = 1 : nCol
        hAx = hGrid{iRow, iCol};

        % Axis limits
        if iRow == 1
            if ~isempty(xLim{iCol})
                hAx.XAxis.Limits = xLim{iCol};
            end
        elseif iRow > 1
            if ~isempty(xLim{iCol})
                hAx.XAxis.Limits = xLim{iCol};
            end
            if ~isempty(yLim{iRow - 1})
                hAx.YAxis.Limits = yLim{iRow - 1};
            end
        end

        % Graphics
        hAx.Box = "on";
        hAx.FontSize = 16;
        hAx.YAxis.Label.FontSize = 20;
        hAx.XAxis.Label.FontSize = 20;

        % Legend
        axPos = get(hAx, 'position');
        hLgd = get(hAx, 'legend');
        if ~isempty(hLgd)
            hLgd.Location = 'south';
            hLgd.Position(2) = axPos(2);
            hLgd.Position(1) = axPos(1);
            hLgd.Position(3) = axPos(3);
        end
    end
end

plot_axSize('hFig', hFig, 'szOnly', false,...
    'axWidth', 1100, 'axHeight', 1100, 'flgPos', false);

% Save
fname = sprintf('MEA~CorrHist(Selective Rows)');
lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat'});











