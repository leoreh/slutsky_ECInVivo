
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
%  PER FILE ANALYSIS
%  ========================================================================

% Files
basepaths = [mcu_basepaths('mea_bac'), mcu_basepaths('mea_mcuko')];
vars = {'mea', 'fr'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);
nFiles = length(basepaths);

% Params
isiVal = 0.015;
binSize = 60;

winExp = [0, 9 * 60]  * 60;
winBsl = [1, 70] * 60;
winSs = [7 * 60, 9 * 60 - 5] * 60;
winTrough = [4 * 60 + 10, 4.5 * 60] * 60;
winCalc = [winBsl; winSs; winTrough];

close all
for iFile = 1 : nFiles

    % File
    basepath = basepaths{iFile};
    [~, basename] = fileparts(basepath);
    cd(basepath);

    % % Organize raw spike times
    % files = dir('*sorted*');
    % mea = mea_orgNex('fname', files.name, 'basepath', basepath, ...
    %     'flgForce', true, 'flgSave', true, 'flgPlot', false);

    % Limit to experimental window
    spktimes = v(iFile).mea.spktimes;
    lastspike = max(cellfun(@max, spktimes, 'uni', true));
    spktimes = cellfun(@(x) x(x >= winExp(1) & x <= winExp(2)), ...
        spktimes, 'UniformOutput', false);

    % % Firing rate
    % fr = mea_frPrep(spktimes, 'binSize', binSize, ...
    %     'flgSave', true, 'flgPlot', true);

    % Burst detection
    brst = brst_detect(spktimes, ...
        'minSpks', 3, ...
        'maxISI_start', isiVal, ...
        'maxISI_end', isiVal * 2, ...
        'minDur', 0.005, ...
        'minIBI', 0.1, ...
        'flgForce', true, 'flgSave', true, 'flgPlot', true);

    % Burst temporal dynamics
    % dyn = brst_dynamics(brst, spktimes, 'binSize', 60, 'kernelSD', 300, ...
    %     'binSize', binSize, 'flgSave', true, 'flgPlot', false);

    % Burst statistics
    stats = brst_stats(brst, spktimes, 'winCalc', winCalc, ...
        'flgSave', true);

    % % FR recovery
    % % fr = v(iFile).fr;
    % rcv = mea_frRecovery(fr.t, fr.fr, ...
    %     'idxTrough', fr.info.idxTrough, ...
    %     'binSize', binSize, 'flgSave', true, 'flgPlot', true);

    % % FR model fit an recovery
    % frFit = mea_frFit(fr.fr, fr.t, 'FilterLen', [], ...
    %     'idxTrough', fr.info.idxTrough, ...
    %     'flgPlot', true, 'flgSave', true);
    % rcvMdl = mea_frRecovery(fr.t, frFit.frMdl, ...
    %     'idxTrough', fr.info.idxTrough, ...
    %     'binSize', binSize, 'flgSave', false, 'flgPlot', false);
    % save(fullfile(basepath, [basename, '.frRcv_mdl.mat']), 'rcvMdl', '-v7.3');

    % Tranfer function spikes to Ca2+
    % ca = spk2ca(spktimes, 'winCalc', [0, Inf], ...
    %     'flgPlot', true, 'flgSave', true);

end

% Dimensionality
% mcu_dim('type', 'vitro')


%% ========================================================================
%  LOAD TABLE
%  ========================================================================

[tbl, xVec, basepaths, v] = mcu_tblMea(basepaths, v);


%% ========================================================================
%  PLOTS
%  ========================================================================


tblGUI_xy(xVec, tbl);

tblGUI_scatHist(tbl, 'xVar', 'bFrac', 'yVar', 'rcvTime', 'grpVar', 'Group');

tblGUI_bar(tbl, 'yVar', 'bFrac', 'xVar', 'Group');

tblGUI_raster(tbl, 'grpVar', 'Name', 'grpVal', 'ctrl1')





%% ========================================================================
%  LME - SPIKING
%  ========================================================================

tblLme = tbl;

% Fit
varRsp = 'fr';
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
%  PLOT CORRELATIONS
%  ========================================================================


% -------------------------------------------------------------------------
% SELECTIVE (updated sep-25)
lData = tbl_transform(lmeData, 'varsInc', {'frBsl', 'BSpks'}, 'flgZ', false,...
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











