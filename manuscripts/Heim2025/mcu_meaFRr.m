%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PRC Params
clear prcParams
prcParams.winLim = [0 70 * 60];        % Analysis window [s]
prcParams.binSize = 0.001;             % 1ms bins
prcParams.gkHw = 0.012;                % 12ms sigma
prcParams.winStpr = 1.0;               % 1s window
prcParams.nShuffles = 1000;            % Number of shuffles
prcParams.spkLim = 2000;
prcParams.shuffleMet = 'raster';

% Files
% basepaths = [mcu_sessions('mea_mk801')];
basepaths = [mcu_sessions('mea_bac'), mcu_sessions('mea_mcuko')];
nFiles = length(basepaths);
vars = {'mea', 'st_metrics'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

% Analysis Params
bslLim = [5 70] * 60;        
ssLim = [7 * 60, 9 * 60 - 5] * 60;
troughLim = [4 * 60 + 10, 4.5 * 60] * 60;
stWin = {bslLim, ssLim, troughLim};
expLim = [0, 9 * 60]  * 60;

% expLim = [0, Inf];

for iFile = 1 : nFiles

    basepath = basepaths{iFile};
    cd(basepath);
    
    % files = dir('*sorted*');
    % mea = mea_orgNex('fname', files.name, 'basepath', pwd, 'forceL', false);
    spktimes = v(iFile).mea.spktimes;
    
    % % --- Firing Rate Recovery
    frr = mea_frr(spktimes, 'winLim', expLim,...
        'flgSave', true, 'flgPlot', false, 'flgForce', false);

    % % --- Spike timing metrics
    % st = spktimes_metrics('spktimes', spktimes, 'sunits', [],...
    %     'bins', stWin, 'flg_force', true, 'flg_save', true, 'flg_all', false);
    % 
    % % --- Bursts
    % brst = spktimes_meaBrst(spktimes, 'binsize', [], 'isiThr', 0.02,...
    %     'minSpks', 2, 'flgSave', true, 'flgForce', true, 'bins', stWin);

    % % --- Population Coupling
    % [prc] = prCoupling(spktimes, prcParams, 'flgSave', true);
    % prCoupling_plot(prc, 'basepath', basepath, 'flgSaveFig', true);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TABLE PREPARATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data from both groups
grps = {'mea_bac'; 'mea_mcuko'};
grpLbls = {'Control'; 'MCU-KO'};
vars = {'frr', 'st_brst'};

% Choose between model-based (mdl) or model-free (mdlF) metrics
mdlPrfx = 'frr.mdl';
mdlPrfx = 'frr.mdlF';

clear varMap
varMap.uGood      = 'frr.uGood';
varMap.frBsl      = [mdlPrfx, '.frBsl'];
varMap.frSs       = [mdlPrfx, '.frSs'];  
varMap.frTrough   = [mdlPrfx, '.frTrough'];  
varMap.pertDepth  = [mdlPrfx, '.pertDepth'];
varMap.uRcv       = [mdlPrfx, '.uRcv'];
varMap.uPert      = [mdlPrfx, '.uPert'];
varMap.rcvTime    = [mdlPrfx, '.rcvTime'];
varMap.bslTime    = [mdlPrfx, '.bslTime'];
varMap.rcvErr     = [mdlPrfx, '.rcvErr'];
varMap.rcvGain    = [mdlPrfx, '.rcvGain'];
varMap.rcvWork    = [mdlPrfx, '.rcvWork'];
varMap.rcvSlope   = [mdlPrfx, '.normSlope'];
varMap.spkDfct    = [mdlPrfx, '.spkDfct'];
varMap.rcvDiff    = [mdlPrfx, '.rcvDiff'];
varMap.rcvFit     = ['frr.frFit.rsquare'];
varMap.BSpks      = 'brst.bspks';
varMap.brBsl      = 'brst.rate';
% varMap.BMiz      = 'st.mizuseki';    % keeping st_metrics removes many rows
% varMap.BRoy      = 'st.royer';
% varMap.prc        = 'prc.prc0_norm';

% Specific overrides
varMap.pertDepth  = 'frr.mdlF.pertDepth';
varMap.uRcv       = 'frr.mdlF.uRcv';
varMap.spkDfct    = 'frr.mdlF.spkDfct';
varMap.rcvGain    = 'frr.mdlF.rcvGain';
varMap.rcvWork    = 'frr.mdl.rcvWork';
varMap.rcvTime    = 'frr.mdl.rcvTime';
varMap.bslTime    = 'frr.mdl.bslTime';

clear tblCell
for iGrp = 1 : length(grps)
    basepaths = mcu_sessions(grps{iGrp});
    
    % Load data for this group
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);
    
    % Prepare tag structures for v2tbl
    tagAll.Group = grpLbls{iGrp};
    tagFiles.Name = get_mname(basepaths);
    
    % Create table using new flexible approach
    tblCell{iGrp} = v2tbl('v', v, 'varMap', varMap,...
        'tagFiles', tagFiles, 'tagAll', tagAll, 'idxCol', 1);
end
tbl = vertcat(tblCell{:});

% Organize for analysis
lmeData = tbl(tbl.uGood, :);
lmeData.uGood = [];
lmeData = rmmissing(lmeData);   

% Convert time to min, specific for mea x3 data reduction
lmeData.rcvTime = lmeData.rcvTime * 3 / 60 / 60;
lmeData.bslTime = lmeData.bslTime * 3 / 60 / 60;


% -------------------------------------------------------------------------

lData = tbl_transform(lmeData, 'varsInc', {'frBsl', 'BSpks'}, 'flgZ', false,...
    'skewThr', 0.1, 'varsGrp', {'Group'}, 'flgLog', true);
clr = mcu_clr();

% Remove units that didn't recover from rcvTime
idxUnits = lmeData.uRcv;
lData.rcvTime(~idxUnits) = nan;
lData.bslTime(~idxUnits) = nan;


varsInc = {'frBsl', 'BSpks', 'pertDepth',...
    'spkDfct', 'rcvWork', 'rcvTime', 'bslTime'};
[hFig, ~, hGrid] = plot_corrHist(lData, 'varsInc', varsInc,...
    'grpIdx', 'Group', 'clrGrp', clr.grp, 'thrOut', 100);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT CORRELATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% -------------------------------------------------------------------------
% SELECTIVE 
% Prepare data
lData = tbl_transform(lmeData, 'varsInc', {'frBsl', 'BSpks'}, 'flgZ', false,...
    'skewThr', 0.1, 'varsGrp', {'Group'}, 'flgLog', true);
clr = mcu_clr();

% Remove units that didn't recover from rcvTime
idxUnits = lmeData.uRcv;
lData.rcvTime(~idxUnits) = nan;

varsInc = {'BSpks', 'spkDfct', 'rcvWork', 'rcvTime'};
clear txtLbls
txtLbls{1} = 'Baseline P(S\inB)';
txtLbls{2} = 'Spk. Deficit';
txtLbls{3} = 'Rcv. Ratio';
txtLbls{4} = 'Rcv. Time (hr)';
lData = lData(:, ['Group', varsInc]);
for iVar = 1 : length(varsInc)
    lData.Properties.VariableNames(2 : end) = txtLbls;
end

varsRow = txtLbls(1);
varsCol = txtLbls(2 : 4);
[hFig, ~, hGrid] = plot_corrHist(lData, 'varsRow', varsRow,...
    'varsCol', varsCol, 'grpIdx', 'Group', 'clrGrp', clr.grp, 'thrOut', 100);
plot_axSize('hFig', hFig, 'szOnly', false,...
    'axWidth', 1200, 'axHeight', 600, 'flgPos', true);

% Update rcvTime scatter
yData = lData.(txtLbls{1});
xData = lData.(txtLbls{4});
[rho, pval] = corr(yData, xData, 'Type', 'Spearman', 'Rows', 'complete');
hAx = hGrid{2, 3};
hSct = scatter(hAx, xData, yData, 10, 'filled', ...
    'MarkerFaceColor', clr.grp(1, 1:3), 'MarkerFaceAlpha', 0.7);   
if pval < 0.0001
    pStr = 'p < 0.0001';
else
    pStr = sprintf('p = %.4f', pval);
end
txtLgd = sprintf('\\rho = %.2f, %s', rho, pStr);
hLgd = legend(hAx, txtLgd, 'Interpreter', 'tex', 'Location', 'south');
hLgd.Box = "on";
xlabel(hAx, txtLbls{4}, 'Interpreter', 'tex', 'FontSize', 20)
ylabel(hAx, '')
yticks(hAx, [])

% Update rcvTime histogram
hAx = hGrid{1, 3}; 
histogram(hAx, xData, 'Normalization', 'pdf', ...
    'FaceColor', clr.grp(1, 1:3), 'FaceAlpha', 0.7, 'EdgeColor', 'none');
hold(hAx, 'on');
meanVal = mean(xData, 'omitnan');
yLimit = ylim(hAx);
plot(hAx, [meanVal meanVal], yLimit, '--', 'Color', 'k', 'LineWidth', 3);
title(hAx, '');
set(hAx, 'XTickLabel', [], 'YTickLabel', []);
set(hAx, 'box', 'off', 'XColor', 'none', 'YColor', 'none');

% Update BSpks labels
hAx = hGrid{2, 1}; 
yticks(hAx, [-2, -1, 0])
yticklabels(hAx, {'0.01', '0.1', '1'});

[nRow, nCol] = size(hGrid);
clear xLim yLim
xLim = cell(nCol - 1, 1);
xLim{1} = [-2, 6];
xLim{2} = [0.5, 1.5];
xLim{3} = [0, 20];
yLim = cell(nRow - 1, 1);
yLim{2} = [-3, 0];

for iRow = 1 : nRow
    for iCol = 1 : nCol
        hAx = hGrid{iRow, iCol};
        
        % Axis limits
        if iCol < nCol
            if ~isempty(axLim{iRow})
                hAx.XAxis.Limits = xLim{iCol};
            end
        end       
        if iRow > 1
             if ~isempty(axLim{iRow})
                hAx.YAxis.Limits = yLim{iRow};
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
    'axWidth', 1200, 'axHeight', 600, 'flgPos', true);

% Save
fname = sprintf('MEA ~ CorrHist(Selective)');
lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat'});





% -------------------------------------------------------------------------
% ALL VARS
% Prepare data
lData = tbl_transform(lmeData, 'varsInc', {'frBsl', 'BSpks'}, 'flgZ', false,...
    'skewThr', 0.1, 'varsGrp', {'Group'}, 'flgLog', true);

varsInc = {'frBsl', 'BSpks', 'pertDepth',...
    'spkDfct', 'rcvWork'};
txtLbls{1} = 'Baseline FR (Hz)';
txtLbls{2} = 'Baseline P(S\inB)';
txtLbls{3} = 'Pert. Depth';
txtLbls{4} = 'Spk. Deficit';
txtLbls{5} = 'Rcv. Ratio';
lData = lData(:, ['Group', varsInc]);
for iVar = 1 : length(varsInc)
    lData.Properties.VariableNames(2 : end) = txtLbls;
end

[hFig, ~, hGrid] = plot_corrHist(lData, 'varsInc', txtLbls,...
    'grpIdx', 'Group', 'clrGrp', clr.grp, 'thrOut', 100);
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'square', 'axHeight', 800);

[nRow, nCol] = size(hGrid);
axLim = cell(nRow, 1);
axLim{1} = [-1, 1];
axLim{2} = [-2, 0];
axLim{3} = [0, 15];
axLim{4} = [-2.5, 5];
axLim{5} = [0.5, 1.5];

% Update BSpks tick labels
gridCoords = [2, 1; 2, 5; 5, 2; 1, 2];
hAx = hGrid{2, 1};
yticks(hAx, [-2, -1, 0])
yticklabels(hAx, {'0.01', '0.1', '1'});
hAx = hGrid{2, 5};
yticks(hAx, [-2, -1, 0])
yticklabels(hAx, {'0.01', '0.1', '1'});
hAx = hGrid{5, 2};
xticks(hAx, [-2, -1, 0])
xticklabels(hAx, {'0.01', '0.1', '1'});

% Update FR tick labels
hAx = hGrid{5, 1};
xticks(hAx, [-1, 0, 1])
xticklabels(hAx, {'0.1', '1', '10'});
hAx = hGrid{1, 5};
yticks(hAx, [-1, 0, 1])
yticklabels(hAx, {'0.1', '1', '10'});

for iRow = 1 : nRow
    for iCol = 1 : nCol
        hAx = hGrid{iRow, iCol};
        if ~isempty(axLim{iRow})
            hAx.XAxis.Limits = axLim{iCol};
        end
        if any(arrayfun(@(h) ...
                isa(h, 'matlab.graphics.chart.primitive.Histogram'), ...
                hAx.Children))

            continue
        end
        if ~isempty(axLim{iRow})
            % Add space for legend
            hAx.YAxis.Limits = [axLim{iRow}(1) - 0.5, axLim{iRow}(2)];
        end
        
        % Graphics
        hAx.Box = "on";
        hAx.FontSize = 16;
        hAx.YAxis.Label.FontSize = 20;
        hAx.XAxis.Label.FontSize = 20;
        if iCol == nCol
            hAx.YAxis.Label.Rotation = 270;
        end

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

plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'square',...
    'axHeight', 1500, 'axWidth', 1500, 'flgPos', false);

% Save
fname = sprintf('MEA ~ CorrHist(All)');
lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LINEAR MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% (*) REMARKS
% Since the primary question concerns genotypes, only test interactions
% between predictors and group (rather then between predictors). Due to
% collinearity, I only use BSpks as a measure of brustiness because it is
% the only one not correlation with bslFr.

% I only apply z score when using contineous predictors for easier
% interpretation. with only grouping variables (when means are compared) I
%  want to keep the original units of the response variables. Same goes for
%  log transform. 

% Insisted to use rcvGain and spkDfct on logarithm scale so they can be
% analyzed with fitlme. rcvTime (and MF) probably still needs glme. 

% BslFr is positively correlated with recovery time (model and model free).
% This makes sense considering most units drop to zero, and thus despite
% normalizing the target value, it is still largely influenced by BslFr. A
% similar result is obtain for SpkDftc. Because of this, and because it
% does not predict uRcv, it is omitted from subsequent models

% For recovery time, the gamma distribution used the reciprocal link which
% means that a positive coefficient decreases the response

% Recovery slope is impossible with raw firing rates and for the model
% fits, it depends too much on the selected model (e.g. sigmoid vs.
% exponential).

% A discripency between model-based and model-free parameters is that in
% the latter there is no correlation between frBsl and pertDepth because
% many values are clamped to c. Hence pertDepth should only be from the
% model. 

% Recovery time includes units that reached their threshold value but
% didn't manage to maintain it, i.e. 

% -------------------------------------------------------------------------
% PREPS

% Recovered units
idxUnits = lmeData.uRcv;
lmeMdl = {};

% List of possible predictors
listPrdct = {'pertDepth', 'frBsl', 'BSpks', 'PRC', 'Group', '(1|Name)', 'brBsl'}; 
listRspns = {'uRcv', 'rcvGain', 'rcvWork', 'rcvErr', 'bslTime', 'spkDfct', 'rcvSlope'}; 

% Z score predictors
zlData = tbl_transform(lmeData, 'varsExc', listRspns, 'flgZ', true,...
    'skewThr', 2, 'varsGrp', {'Group'}, 'flgLog', true);

% -------------------------------------------------------------------------
% RECOVERY PROBABILITY
iPr = [1, 2, 3, 5, 6];
frml = sprintf('%s ~ %s', listRspns{1}, strjoin(listPrdct([iPr]), ' + '));
lmeMdl{end + 1} = fitglme(zlData, frml, 'FitMethod', 'REMPL',...
    'Distribution', 'Binomial');

% -------------------------------------------------------------------------
% RECOVERY WORK 
frml = [listRspns{3}, ' ~ frBsl + BSpks * Group + (1|Name)'];
lmeMdl{end + 1} = fitlme(zlData, frml, 'FitMethod', 'REML');

% -------------------------------------------------------------------------
% RECOVERY TIME 
% only units who recovered, from both groups combined
frml = [listRspns{5}, ' ~ frBsl + BSpks + (1|Name)'];
lmeMdl{end + 1} = fitglme(zlData(idxUnits, :), frml, 'FitMethod', 'REMPL',...
    'Distribution', 'Gamma');

% -------------------------------------------------------------------------
% SPIKE DEFICIT 
frml = [listRspns{6}, ' ~ pertDepth + frBsl + BSpks * Group + (1|Name)'];
lmeMdl{end + 1} = fitlme(zlData, frml, 'FitMethod', 'REML');

% % -------------------------------------------------------------------------
% % RECOVERY ERROR 
% % only units who recovered, from both groups combined
% frml = [listRspns{4}, ' ~ pertDepth + BSpks + frBsl + Group + (1|Name)'];
% lmeMdl{end + 1} = fitglme(zlData(idxUnits, :), frml, 'FitMethod', 'REMPL',...
%     'Distribution', 'Gamma');

% -------------------------------------------------------------------------
% SAVE MODELS
fname = 'MEA ~ frrMdl';
for iMdl = 1 : length(lmeMdl)
    sheetNames = {'LME_Data', lmeMdl{iMdl}.ResponseName};
    lme_save('fname', fname, 'frmt', {'mat', 'xlsx'},...
        'lmeData', lmeData, 'lmeMdl', lmeMdl{iMdl}, 'sheetNames', sheetNames)
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FR DISTRIBUTION DURING RECOVERY 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepaths = mcu_sessions('mea_bac');
vars = {'frr', 'st_brst'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

frr = catfields([v(:).frr], 1);
fr = frr.fr;
uGood = frr.uGood;
idxTrough = mode(frr.mdlF.idxTrough);

frr = catfields([v(:).frr], 2);
t = (frr.t(:, 1) - frr.t(idxTrough, 1)) * 3 / 60 / 60;

win = [-4, -3; 5, 6; 8, 9];
clear xData
close all
plot_axSize('szOnly', false, 'axShape', 'square',...
    'axHeight', 400, 'axWidth', 400, 'flgPos', true);
hold on
for iWin = 1 : size(win, 1)
    [~, winIdx(1)] = min(abs(t - win(iWin, 1)));
    [~, winIdx(2)] = min(abs(t - win(iWin, 2)));

    xData = log10(mean(fr(uGood, winIdx), 2, 'omitnan'));
    xData = rmoutliers(xData, 'percentiles', [2, 98]);
    nBins = 20;
    hHst = histogram((xData), nBins, 'Normalization', 'probability');
    hHst.EdgeColor = 'none';
    hHst.FaceAlpha = 0.5;
end
legend


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLOPE VS BURSTINESS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except clr
basepaths = mcu_sessions('mea_bac');
vars = {'frr', 'st_brst'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

frr = catfields([v(:).frr], 1);
brst = catfields([v(:).brst], 2);
uGood = frr.uGood;
idxTrough = mode(frr.frFit.idxTrough);
bslWin = 1 : idxTrough - 5;

% Calculate robust slope for each uGood neuron
slopes = nan(size(uGood));
for iUnit = 1:length(uGood)
    if uGood(iUnit)
        t = 1:length(frr.fr(iUnit, bslWin));
        [b, ~] = robustfit(t, frr.fr(iUnit, bslWin));
        slopes(iUnit) = b(2);
    end
end
bslDrift = slopes ./ frr.mdlF.frBsl * length(bslWin);

% Plot correlation
figure;
scatter(log10(brst.bspks(uGood)), bslDrift(uGood));
xlabel('Normalized Slope');
ylabel('Spike Percentage in Bursts');
title('Correlation between Normalized Slope and Burstiness');
grid on;

