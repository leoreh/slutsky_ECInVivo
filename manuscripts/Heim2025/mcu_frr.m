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
    % frr = mea_frr(spktimes, 'winLim', expLim,...
    %     'flgSave', true, 'flgPlot', false, 'flgForce', false);

    % --- Spike timing metrics
    st = spktimes_metrics('spktimes', spktimes, 'sunits', [],...
        'bins', stWin, 'flg_force', true, 'flg_save', true, 'flg_all', false);

    % --- Bursts
    brst = spktimes_meaBrst(spktimes, 'binsize', [], 'isiThr', 0.02,...
        'minSpks', 2, 'flg_save', true, 'flg_force', true, 'bins', stWin);

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
vars = {'frr', 'st_metrics', 'st_brst'};
vars = {'frr', 'st_brst'};

% Choose between model-based (mdl) or model-free (mdlF) metrics
mdlPrfx = 'frr.mdl';
mdlPrfx = 'frr.mdlF';

clear varMap
varMap.uGood      = 'frr.uGood';
varMap.frBsl      = [mdlPrfx, '.frBsl'];
varMap.frRcv      = [mdlPrfx, '.frSs'];  
varMap.pertDepth  = [mdlPrfx, '.pertDepth'];
varMap.uRcv       = [mdlPrfx, '.uRcv'];
varMap.rcvTime    = [mdlPrfx, '.rcvTime'];
varMap.rcvErr     = [mdlPrfx, '.rcvErr'];
varMap.rcvGain    = [mdlPrfx, '.rcvGain'];
varMap.rcvSlope   = [mdlPrfx, '.normSlope'];
varMap.spkDfct    = [mdlPrfx, '.spkDfct'];
varMap.rcvDiff    = [mdlPrfx, '.rcvDiff'];
varMap.brPct      = 'brst.spkprct';
varMap.brBsl      = 'brst.rate';
% varMap.brMiz      = 'st.mizuseki'; 
% varMap.brRoy      = 'st.royer';
% varMap.prc        = 'prc.prc0_norm';

% Specific overrides
varMap.pertDepth  = 'frr.mdl.pertDepth';
varMap.uRcv       = 'frr.mdl.uRcv';
varMap.spkDfct    = 'frr.mdlF.spkDfct';
varMap.rcvTime    = 'frr.mdl.rcvTime';

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
lmeData = rmmissing(lmeData);   % Note, keeping st_metrics removes many rows

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT CORRELATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear clrGrp
clrGrp(1, :) = [0.3 0.3 0.3, 0.4];
clrGrp(2, :) = [0.784 0.667 0.392, 0.4];

% Log transform for plotting
varsInc = {'frBsl', 'brPct', 'rcvErr', 'pertDepth', 'brBsl'};
lData = tbl_transform(lmeData, 'varsInc', varsInc, 'flgZ', false,...
    'skewThr', 0.1, 'varsGrp', {'Group'}, 'flgLog', true);

% Remove units that didn't recover from rcvTime
idxUnits = lmeData.uRcv;
lData.rcvTime(~idxUnits) = nan;
lData.rcvErr(~idxUnits) = nan;
lData.rcvSlope(~idxUnits) = nan;

% Plot correlations
varsCol = {'spkDfct', 'rcvGain', 'rcvTime'};
varsRow = {'brPct'};
[hFig, ~, hGrid] = plot_corrHist(lData, 'varsRow', varsRow, 'varsCol', varsCol,...
    'grpIdx', 'Group', 'clrGrp', clrGrp, 'thrOut', 95);
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'wide', 'axHeight', 600);

% Plot correlations
varsRow = {'spkDfct', 'rcvGain'};
varsCol = {'brPct', 'frBsl'};
[hFig, ~, hGrid] = plot_corrHist(lData, 'varsRow', varsRow, 'varsCol', varsCol,...
    'grpIdx', 'Group', 'clrGrp', clrGrp, 'thrOut', 100);
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'wide', 'axHeight', 600);

hGrid{1, 1}

varsInc = {'frBsl', 'brPct', 'pertDepth', 'spkDfct', 'rcvGain', 'rcvTime'};
[hFig, ~, hGrid] = plot_corrHist(lData, 'varsInc', varsInc,...
    'grpIdx', 'Group', 'clrGrp', clrGrp, 'thrOut', 100);

plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'wide', 'axHeight', 700);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LINEAR MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% (*) REMARKS
% Since the primary question concerns genotypes, only test interactions
% between predictors and group (rather then between predictors). Due to
% collinearity, I only use brPct as a measure of brustiness because it is
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
listPrdct = {'pertDepth', 'frBsl', 'brPct', 'PRC', 'Group', '(1|Name)', 'brBsl'}; 
listRspns = {'uRcv', 'rcvGain', 'rcvErr', 'rcvTime', 'spkDfct', 'rcvSlope'}; 

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
% RECOVERY GAIN 
% No point including pertDepth because it's part of definition of Gain
frml = [listRspns{2}, ' ~ frBsl + brPct * Group + (1|Name)'];
lmeMdl{end + 1} = fitlme(zlData, frml, 'FitMethod', 'REML');

% -------------------------------------------------------------------------
% RECOVERY TIME 
% only units who recovered, from both groups combined
frml = [listRspns{4}, ' ~ pertDepth + frBsl + brPct + (1|Name)'];
% lmeMdl{end + 1} = fitglme(zlData(idxUnits, :), frml, 'FitMethod', 'REMPL',...
%     'Distribution', 'Gamma');
lmeMdl{end + 1} = fitlme(zlData, frml, 'FitMethod', 'REML');

% -------------------------------------------------------------------------
% SPIKE DEFICIT 
frml = [listRspns{5}, ' ~ pertDepth + frBsl + brPct * Group + (1|Name)'];
lmeMdl{end + 1} = fitlme(zlData, frml, 'FitMethod', 'REML');

% % -------------------------------------------------------------------------
% % RECOVERY ERROR 
% % only units who recovered, from both groups combined
% frml = [listRspns{3}, ' ~ pertDepth + brPct + frBsl + Group + (1|Name)'];
% lmeMdl{end + 1} = fitglme(zlData(idxUnits, :), frml, 'FitMethod', 'REMPL',...
%     'Distribution', 'Gamma');
% 
% frml = [listRspns{6}, ' ~ pertDepth + brPct + frBsl + Group + (1|Name)'];
% lmeMdl{end + 1} = fitlme(zlData, frml, 'FitMethod', 'REML');

% -------------------------------------------------------------------------
% SAVE MODELS
fname = 'MEA ~ frrMdl';
for iMdl = 1 : length(lmeMdl)
    sheetNames = {'LME_Data', lmeMdl{iMdl}.ResponseName};
    lme_save('fname', fname, 'frmt', {'mat', 'xlsx'},...
        'lmeData', lmeData, 'lmeMdl', lmeMdl{iMdl}, 'sheetNames', sheetNames)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARE BASELINE VS RECOVERY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Load data from both groups
grps = {'mea_bac'; 'mea_mcuko'};
grpLbls = {'Control'; 'MCU-KO'};
vars = {'frr', 'st_metrics', 'st_brst'};
vars = {'frr', 'st_brst'};

% Create varMap based on time point
clear varMap
varMap.uGood      = 'frr.uGood';
% varMap.BrMiz      = 'st.mizuseki';
% varMap.BrRoy      = 'st.royer';
varMap.BrPct      = 'brst.spkprct';

% Override the source prefix (in main analysis)
mdlPrfx = 'frr.mdlF';

cnt = 1;
clear tblCell
for iGrp = 1 : length(grps)
    basepaths = mcu_sessions(grps{iGrp});
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);

    % Create consistent UnitIDs for this group (reset between groups). Use
    % frr for grabbing the number of units.
    frr = catfields([v(:).frr], 1);
    [nUnits, ~] = size(frr.fr);
    unitIDs = (1:nUnits)';

    for iCol = 1 : 3

        % Fr maps to different sources based on time point
        switch iCol
            case 1
                varMap.Fr  = [mdlPrfx, '.frBsl'];
                tagAll.Time = 'BSL';
            case 2
                varMap.Fr  = [mdlPrfx, '.frTrough'];
                tagAll.Time = 'BAC 1h';

            case 3
                varMap.Fr  = [mdlPrfx, '.frSs'];
                tagAll.Time = 'BAC 24h';
        end

        % Prepare tag structures for v2tbl
        tagAll.Group = grpLbls{iGrp};
        tagFiles.Name = get_mname(basepaths);

        tempTbl = v2tbl('v', v, 'varMap', varMap, 'tagFiles',...
            tagFiles, 'tagAll', tagAll, 'idxCol', iCol);
        
        % Override the UnitID to ensure consistency across time points
        % Add group offset to make UnitIDs unique between groups
        groupOffset = (iGrp - 1) * 10000;  % 10000 units per group
        tempTbl.UnitID = groupOffset + unitIDs;

        tblCell{cnt} = tempTbl;
        cnt = cnt + 1;
    end
end
tbl = vertcat(tblCell{:});
tbl = sortrows(tbl, 'Group');


% -------------------------------------------------------------------------
% Bursty

% Organize for analysis
lmeData = tbl(tbl.uGood, :);
lmeData.uGood = [];
lmeData = rmmissing(lmeData);

% Normalize to baseline
nData = tbl_transform(lmeData, 'varsInc', {varRsp}, 'flgZ', false,...
    'skewThr', 2, 'varsGrp', {'Group'}, 'varNorm', 'Time',...
    'flgLog', false, 'flgNorm', false);
nData.BrPct = nData.BrPct + min(nData.BrPct(nData.BrPct > 0)) / 2;

% Remove unused categorical levels
% nData(nData.Time == 'BaselBACine', :) = [];
nData(nData.Time == 'BAC 1h', :) = [];
if ismember('Time', nData.Properties.VariableNames)
    actualTimeLevels = unique(nData.Time);
    nData.Time = categorical(nData.Time, actualTimeLevels);
end
grpstats(lmeData, {'Group', 'Time'});

% run lme
varRsp = 'BrPct';
frml = [varRsp, ' ~ Group * Time + (1|Name)'];
lmeCfg.contrasts = [1 : 7];
lmeCfg.distribution = 'Gamma';
[lmeStats, lmeMdl] = lme_analyse(nData, frml, lmeCfg);

% plot
hFig = lme_plot(nData, lmeMdl, 'lmeStats', lmeStats,...
    'ptype', 'bar', 'axShape', 'square', 'idxRow', [4, 3]);

% Update labels
hAx = gca;
ylabel(hAx, 'Burstiness P(Spk\inBrst)', 'Interpreter', 'tex')
xlabel(hAx, '')
title(hAx, '')
hAx.Legend.Location = 'southeast';

plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'square',...
    'axHeight', 300)



% -------------------------------------------------------------------------
% Firing Rate

% lmeData(lmeData.Time == 'Acute', :) = [];

% Organize for analysis
lmeData = tbl(tbl.uGood, :);
lmeData.uGood = [];
lmeData = rmmissing(lmeData);

% Normalize to baseline
nData = tbl_transform(lmeData, 'varsInc', {varRsp}, 'flgZ', false,...
    'skewThr', 2, 'varsGrp', {'Group'}, 'varNorm', 'Time',...
    'flgLog', false, 'flgNorm', true);
nData.Fr = nData.Fr + min(nData.Fr(nData.Fr > 0)) / 2;

% Remove unused categorical levels
nData(nData.Time == 'Baseline', :) = [];
nData(nData.Time == 'Acute', :) = [];
if ismember('Time', nData.Properties.VariableNames)
    actualTimeLevels = unique(nData.Time);
    nData.Time = categorical(nData.Time, actualTimeLevels);
end
grpstats(lmeData, {'Group', 'Time'});

% run lme
varRsp = 'Fr';
frml = [varRsp, ' ~ Group + (1|Name)'];
lmeCfg.contrasts = [];
lmeCfg.distribution = 'Gamma';
[lmeStats, lmeMdl] = lme_analyse(nData, frml, lmeCfg);

% plot
hFig = lme_plot(nData, lmeMdl, 'lmeStats', lmeStats,...
    'ptype', 'bar', 'axShape', 'square', 'idxRow', 2);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FR FIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Load data from both groups
grps = {'mea_bac'; 'mea_mcuko'};
grpLbls = {'Control'; 'MCU-KO'};
vars = {'frr', 'st_brst'};
clear v
for iGrp = 1 : 2
    basepaths = mcu_sessions(grps{iGrp});
    v{iGrp} = basepaths2vars('basepath', basepaths, 'vars', vars);
end

% Get perturbation onset from the first session (should be same for all)
% and create time vector
idxPert = v{iGrp}(1).frr.info.idxPert(1);
binSize = v{iGrp}(1).frr.info.binSize(1);
t = ((1:size(frr.fr, 2)) - idxPert) * binSize / 60 / 60 * 3;
pertTime = t(idxPert);

% Set graphics
[hFig, hAx] = plot_axSize('szOnly', false, 'axShape', 'wide', 'axHeight', 300);
clrGrp(1, :) = [0.3 0.3 0.3];
clrGrp(2, :) = [0.784 0.667 0.392];
hold on;

% Select group to plot (1 = Control, 2 = MCU-KO)
for iGrp = 1 : 2

    frr = catfields([v{iGrp}(:).frr], 1);
    uGood = frr.uGood;
    
    % Select units by burstiness
    brst = catfields([v{iGrp}(:).brst], 2);
    brstData = brst.spkprct(1, uGood)';
    thrBrst = 80;
    uBrst = brstData > prctile(brstData, thrBrst);
    uGood = uGood(uBrst);

    % Limit by FR
    frBsl = frr.mdlF.frBsl;
    thrFr = prctile(frBsl, 80);
    uFr = frBsl > thrFr;
    uGood = uGood & uFr;

    % Calculate MFR from good units
    fr = frr.fr(uGood, :);
    frBsl = frr.mdlF.frBsl(uGood);
    mfr = mean((fr ./ frBsl), 1, 'omitnan')  * 100;        

    % Plot raw MFR
    hPtl = plot_stdShade('dataMat', (fr ./ frBsl) * 100,...
        'xVal', t, 'axh', hAx, ...
        'clr', clrGrp(iGrp, :) + 0.1, 'alpha', 0.5);
    set(hPtl, 'DisplayName', sprintf('%s (n=%d)', grpLbls{iGrp}, sum(uGood)));

end

% Formatting
xlabel('Time (min)');
ylabel('MFR (% BSL)');
legend('Location', 'southeast');
xlim([-4, 20]);
xticks([0 : 6 : 24])
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'wide', 'axHeight', 300);

% Mark perturbation
clrPert = [0.5, 0, 0.5, 0.5];
yLimit = ylim;
yText = yLimit(2);
yLine = yLimit(2) - 0.05 * (yLimit(2) - yLimit(1)); % Line 5% below text
plot([pertTime, 24], [yLine, yLine], 'Color', clrPert,...
    'LineWidth', 4, 'HandleVisibility', 'off');
text(pertTime + 0.5, yText, 'Baclofen', 'FontName', 'Arial', 'FontSize', 16, ...
    'Color', clrPert, 'HorizontalAlignment', 'left');

fname = lme_frml2char('MEA ~ fitMFR', 'rmRnd', true);
lme_save('fh', hFig, 'fname', fname, 'frmt', {'svg', 'mat'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT SAMPLE NEURONS WITH FITS - CONTROL ONLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Figure 2: Sample neurons with fits

% Load data for control group only
grps = {'mea_bac'};
grpLbls = {'Control'};
vars = {'frr'};

% Load control data
basepaths = mcu_sessions(grps{1});
v = basepaths2vars('basepath', basepaths, 'vars', vars);

% Combine all frr structures from control group
frr = catfields([v(:).frr], 1, true);
uGood = frr.mdl.uRcv;
mdlName = frr.mdlName{1};

% Get perturbation onset and bin size from the first session
idxPert = frr.info.idxPert(1);
binSize = frr.info.binSize(1);

% Create time vector in minutes
t = ((1:size(frr.fr, 2)) - idxPert) * binSize / 60 / 60 * 3; 
pertTime = t(idxPert);

% Plot a sample neuron
idxSmpl = [18, 20];
uIdx = 28;

[hFig, hAx] = plot_axSize('szOnly', false, 'axShape', 'wide', 'axHeight', 300);
clrUnits = bone(nSmp);
hold on;

% Plot raw firing rate
fr = frr.fr(uIdx, :);
plot(t, fr, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5,...
    'HandleVisibility', 'off');

% Plot fit curve (only recovery portion) with legend entry
fitCurve = frr.frMdl(uIdx, :);
[~, idxRcv] = min(fitCurve);
txtMdl = [upper(mdlName{uIdx}(1)), mdlName{uIdx}(2 : end)];
plot(t(idxRcv:end), fitCurve(idxRcv:end), '-', ...
    'Color', 'k', 'LineWidth', 2.5, ...
    'DisplayName', sprintf('Mdl. %s', txtMdl));

% Add asterisk at recovery time
idxRcv = round(frr.mdl.rcvTime(uIdx) / binSize) + frr.mdl.idxRcv(uIdx);
plot(t(idxRcv), fr(idxRcv), 'o', 'MarkerEdgeColor', 'none',...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 12, 'LineWidth', 2, ...
    'DisplayName', sprintf('Rcv. Time: %d min',...
    round(frr.mdl.rcvTime(uIdx) * 3 / 60)));

% Add asterisk at recovery time (model-free)
idxRcv = round(frr.mdlF.rcvTime(uIdx) / binSize) + frr.mdlF.idxRcv(uIdx);
plot(t(idxRcv), fitCurve(idxRcv), 'diamond', 'MarkerEdgeColor', 'none',...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 12, 'LineWidth', 2, ...
    'DisplayName', sprintf('Rcv. Time (mdl): %d min',...
    round(frr.mdlF.rcvTime(uIdx) * 3 / 60)));

% Fill area of spk deficit
tRcv = t(idxPert:end);
frRcv = fr(idxPert:end);
bslLine = ones(size(tRcv)) * frr.mdl.frBsl(uIdx);
tRcv = tRcv(:)';
frRcv = frRcv(:)';
bslLine = bslLine(:)';
fill(hAx, [tRcv, fliplr(tRcv)], [bslLine, fliplr(frRcv)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none',...
    'DisplayName', sprintf('Spk. Deficit: %.2f', frr.mdl.spkDfct(uIdx)));

% Create dummy plot for spike deficit legend entry
% plot(NaN, NaN, 's', 'MarkerFaceColor', [1, 1, 1], ...
%     'MarkerEdgeColor', 'none', 'MarkerSize', 10, ...
%     'DisplayName', sprintf('Rcv Gain: %.2f', frr.mdl.rcvGain(uIdx)));

% Formatting
xlabel('Time (min)');
ylabel('Firing Rate (Hz)');
xlim([-4, 20]);
xticks([0 : 6 : 24]);
% ylim([0 15])
hLgd = legend;
hLgd.Location = "northwest";
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'wide', 'axHeight', 300);

% Mark perturbation
clrPert = [0.5, 0, 0.5, 0.5];
yLimit = ylim;
yText = 95;
yLine = yText - 0.05 * (yLimit(2) - yLimit(1)); % Line 5% below text
plot([pertTime, 24], [yLine, yLine], 'Color', clrPert,...
    'LineWidth', 4, 'HandleVisibility', 'off');
text(pertTime, yText, 'Baclofen', 'FontName', 'Arial', 'FontSize', 16, ...
    'Color', clrPert, 'HorizontalAlignment', 'left');

fname = lme_frml2char('MEA ~ fitUnit', 'rmRnd', true,...
    'sfx', ['_', num2str(uIdx)]);
lme_save('fh', hFig, 'fname', fname, 'frmt', {'svg', 'mat'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: GET_FLDSZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




