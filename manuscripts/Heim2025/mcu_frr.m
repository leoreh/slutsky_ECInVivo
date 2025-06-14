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
basepaths = [mcu_sessions('mea_bac')];
basepaths = [mcu_sessions('mea_mcuko')];
basepaths = [mcu_sessions('mea_mk801')];

nFiles = length(basepaths);
vars = {'mea', 'st_metrics', 'fr'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

% Analysis Params
winLim = [5 * 60 60 * 60];        
expLim = [0, 8 * 60 * 60];
expLim = [0, Inf];

for iFile = 1 : nFiles

    basepath = basepaths{iFile};
    cd(basepath);
    
    % files = dir('*sorted*');
    % mea = mea_orgNex('fname', files.name, 'basepath', pwd, 'forceL', false);
    spktimes = v(iFile).mea.spktimes;
    
    % % --- Firing Rate Recovery
    frr = mea_frr(spktimes, 'winLim', expLim,...
        'flgSave', true, 'flgPlot', true);

    % % --- Spike timing metrics
    % st = spktimes_metrics('spktimes', spktimes, 'sunits', [],...
    %     'bins', {[0, 70 * 60], [6 * 60 * 60, 8 * 60 * 60]}, 'flg_force', true, 'flg_save', true, 'flg_all', false);
    % 
    % % --- Bursts
    % brst = spktimes_meaBrst(spktimes, 'binsize', [], 'isiThr', 0.02,...
    %     'minSpks', 2, 'flg_save', true, 'flg_force', true, 'bins', {[0, 70 * 60], [6 * 60 * 60, 8 * 60 * 60]});

    % % --- Population Coupling
    % [prc] = prCoupling(spktimes, prcParams, 'flgSave', true);
    % prCoupling_plot(prc, 'basepath', basepath, 'flgSaveFig', true);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TABLE PREPARATION - Control vs MCU-KO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Load data from both groups
grps = {'mea_bac'; 'mea_mcuko'};
grpLbls = {'Control'; 'MCU-KO'};
vars = {'frr', 'st_metrics', 'st_brst'};
clear varMap
varMap.uGood      = 'frr.info.uGood';
varMap.FrBsl      = 'frr.frBsl';
varMap.FrRcv      = 'frr.frRcv';
varMap.BrMiz      = 'st.mizuseki'; 
varMap.BrRoy      = 'st.royer';
varMap.BrPct      = 'brst.spkprct';
varMap.PertDepth  = 'frr.pertDepth';
varMap.uRcv       = 'frr.uRcv';
varMap.RcvTime    = 'frr.rcvTime';
varMap.RcvErr     = 'frr.rcvErr';
varMap.RcvChng    = 'frr.rcvChange';
varMap.RcvSlope   = 'frr.normSlope';
varMap.SpkDfct    = 'frr.mf.spkDeficit';
varMap.mfTime     = 'frr.mf.rcvTime';

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT CORRELATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Corr plot. Change to gplotmatrix
varsInc = {'FrBsl', 'BrPct', 'SpkDfct', 'RcvChng', 'RcvTime', 'mfTime'};
grpVec = double(lmeData.Group == 'MCU-KO') + 1;

figure;
grpIdx = lmeData.Group == 'Control';
corrplot(gca, lmeData(grpIdx, varsInc), 'type', 'Spearman',...
    'testR', 'on');
figure;
grpIdx = lmeData.Group == 'MCU-KO';
corrplot(gca, lmeData(grpIdx, varsInc), 'type', 'Spearman',...
    'testR', 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LINEAR MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% (*) REMARKS
% Since the primary question concerns genotypes, only test interactions
% between predictors and group (rather then between predictors). Due to
% collinearity, I only use BrPct as a measure of brustiness because it is
% the only one not correlation with bslFr.

% I only apply z score when using contineous predictors for easier
% interpretation. with only grouping variables (when means are compared) I
%  want to keep the original units of the response variables. Same goes for
%  log transform. 

% Insisted to use rcvChng and SpkDfct on logarithm scale so they can be
% analyzed with fitlme. RcvTime (and MF) probably still needs glme. 

% BslFr is positively correlated with recovery time (model and model free).
% This makes sense considering most units drop to zero, and thus despite
% normalizing the target value, it is still largely influenced by BslFr. A
% similar result is obtain for SpkDftc. Because of this, and because it
% does not predict uRcv, it is omitted from subsequent models

% -------------------------------------------------------------------------
% PREPS

% List of possible predictors
listPrdct = {'PertDepth', 'FrBsl', 'BrPct', 'PRC', 'Group', '(1|Name)'}; 
listRspns = {'RcvErr', 'RcvTime', 'RcvSlope', 'SpkDfct', 'mfTime', 'uRcv', 'RcvChng'}; 

% Indices of units that recovered
idxUnits = lmeData.uRcv;

% Z score predictors
zlData = tbl_transform(lmeData, 'varsExc', listRspns, 'flgZ', true,...
    'skewThr', 2, 'varsGrp', {'Group'}, 'flgLog', true);

% -------------------------------------------------------------------------
% (1) RECOVERY PROBABILITY
iPr = [1, 2, 3, 5, 6];
frml = sprintf('%s ~ %s', listRspns{6}, strjoin(listPrdct([iPr]), ' + '));
mdl = fitglme(zlData, frml, 'FitMethod', 'REMPL',...
    'Distribution', 'Binomial');

fname = lme_frml2char(frml, 'rmRnd', true);
lme_save('fh', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeData', lmeData, 'lmeMdl', mdl)

lmeResults = lme_mdl2tbls(mdl);

% -------------------------------------------------------------------------
% (2) RECOVERY FIDELITY 
% With interaction on burstiness
% For both SpkDfct and RcvErr, larger values = worse
% recovery. Hence a negative estimate indicates better recovery. Note
% SpkDfct correlated with rcvTime (makes sense)
frml = [listRspns{7}, ' ~ PertDepth + BrPct * Group + (1|Name)'];
mdl = fitlme(zlData, frml, 'FitMethod', 'REML');

% -------------------------------------------------------------------------
% (3) RECOVERY SPEED 
% only units who recovered, from both groups combined
frml = [listRspns{5}, ' ~ PertDepth + BrPct + (1|Name)'];
mdl = fitglme(zlData(idxUnits, :), frml, 'FitMethod', 'REMPL',...
    'Distribution', 'Gamma');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARE BASELINE VS RECOVERY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Load data from both groups
grps = {'mea_bac'; 'mea_mcuko'};
grpLbls = {'Control'; 'MCU-KO'};
vars = {'frr', 'st_metrics', 'st_brst'};

cnt = 1; clear tblCell
for iCol = 1 : 2
    for iGrp = 1 : length(grps)
        basepaths = mcu_sessions(grps{iGrp});

        v = basepaths2vars('basepaths', basepaths, 'vars', vars);
        
        % Create varMap based on time point
        clear varMap
        varMap.uGood      = 'frr.info.uGood';
        varMap.BrMiz      = 'st.mizuseki'; 
        varMap.BrRoy      = 'st.royer';
        varMap.BrPct      = 'brst.spkprct';
        
        % FrBsl maps to different sources based on time point
        if iCol == 1
            varMap.FrBsl  = 'frr.mf.frBsl';  % Before: baseline firing rate
        else
            varMap.FrBsl  = 'frr.mf.frSs';   % After: steady-state firing rate
        end
        
        % Prepare tag structures for v2tbl
        tagAll.Group = grpLbls{iGrp};
        tagFiles.Name = get_mname(basepaths);
        
        % Add Time column through tagAll based on time point
        if iCol == 1
            tagAll.Time = 'Before';
        else
            tagAll.Time = 'After';
        end
        
        tempTbl = v2tbl('v', v, 'varMap', varMap, 'tagFiles',...
            tagFiles, 'tagAll', tagAll, 'idxCol', iCol);
        
        tblCell{cnt} = tempTbl;
        cnt = cnt + 1;
    end
end
tbl = vertcat(tblCell{:});

% Organize for analysis
lmeData = tbl(tbl.uGood, :);
lmeData.uGood = [];
lmeData = rmmissing(lmeData);

% FR per unit, WT vs MCU for RS vs FS 
frml = 'FrBsl ~ Group * Time + (1|Name)';
frml = 'BrPct ~ Group * Time + (1|Name)';

mdl = fitglme(lmeData, frml, 'FitMethod', 'REMPL',...
    'Distribution', 'Gamma');

% run lme
frml = 'FrBsl ~ Group * Time + (1|Name)';
lmeCfg.contrasts = 'all';
lmeCfg.distribution = 'gamma';
[lmeStats, lmeMdl] = lme_analyse(lmeData, frml, lmeCfg);

% plot
hFig = lme_plot(lmeData, lmeMdl, 'ptype', 'bar', 'axShape', 'square');







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FR FIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Load data from both groups
grps = {'mea_bac'; 'mea_mcuko'};
grpLbls = {'Control'; 'MCU-KO'};
vars = {'frr'};

% -------------------------------------------------------------------------
% Figure 1: MFR with fit
[hFig, hAx] = plot_axSize('szOnly', false, 'axShape', 'wide', 'axHeight', 300);
clrGrp(1, :) = [0.3 0.3 0.3];
clrGrp(2, :) = [0.784 0.667 0.392];
hold on;

% Select group to plot (1 = Control, 2 = MCU-KO)
for iGrp = 1 : 2
    basepaths = mcu_sessions(grps{iGrp});
    v = basepaths2vars('basepath', basepaths, 'vars', vars);

    % Combine all frr structures from this group
    frr = catfields([v(:).frr], 1);
    uGood = frr.info.uGood;

    % Calculate MFR from good units
    fr = frr.fr(uGood, :);
    mfr = mean((fr ./ frr.frBsl(uGood)), 1, 'omitnan');

    % Get perturbation onset from the first session (should be same for all)
    idxPert = frr.info.pertOnset(1);
    binSize = frr.info.binSize(1);

    % Fit MFR using mea_frFit
    mfrFit = mea_frFit(mfr, idxPert, 'flgPlot', false, 'binSize', binSize);
    [~, idxRcv] = min(mfrFit.fitCurve);

    % Create time vector in minutes
    t = (1:length(mfr)) * binSize / 60 / 60 * 3;  % Convert to minutes
    pertTime = t(idxPert);

    % Plot raw MFR
    hPtl = plot_stdShade('dataMat', (fr ./ frr.frBsl(uGood)), 'xVal', t, 'axh', hAx, ...
        'clr', clrGrp(iGrp, :) + 0.1, 'alpha', 0.5);
    set(hPtl, 'DisplayName', sprintf('%s (n=%d)', grpLbls{iGrp}, sum(uGood)));
    plot(t(idxRcv : end), mfrFit.fitCurve(idxRcv : end), 'b--',...
        'LineWidth', 3, 'Color', clrGrp(iGrp, :) - 0.2, 'HandleVisibility', 'off');

end

% Mark perturbation onset
% Mark perturbation onset
xline(pertTime, 'r--', {'Pert. Onset'}, 'LineWidth', 1,...
    'LabelOrientation', 'horizontal', 'HandleVisibility', 'off', ...
    'LabelHorizontalAlignment', 'left');

% Formatting
xlabel('Time (min)');
ylabel('MFR (Hz)');
legend('show', 'Location', 'northeast');
xlim([0, 24]);
xticks([0 : 6 : 24])
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'wide', 'axHeight', 300);


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
frr = catfields([v(:).frr], 1);
uGood = frr.uRcv;

% Get perturbation onset and bin size from the first session
idxPert = frr.info.pertOnset(1);
binSize = frr.info.binSize(1);

% Create time vector in minutes
t = (1:size(frr.fr, 2)) * binSize / 60 / 60 * 3;  % Convert to minutes
pertTime = t(idxPert);

% Select a sample of good neurons (e.g., first 10 or random selection)
nSmp = min(5, sum(uGood));
idxSmpl = find(uGood);
idxSmpl = idxSmpl(randperm(length(idxSmpl), nSmp));
idxSmpl = [67, 291, 574];

[hFig, hAx] = plot_axSize('szOnly', false, 'axShape', 'wide', 'axHeight', 300);
clrUnits = bone(nSmp);
hold on;

% Plot each sample neuron
for iUnit = 1:length(idxSmpl)
    uIdx = idxSmpl(iUnit);
        
    % Plot raw firing rate
    fr = frr.fr(uIdx, :);
    plot(t, fr, 'Color', [clrUnits(iUnit, :), 0.5], 'LineWidth', 0.5);
    
    % Add asterisk at recovery time
    idxRcv = round(frr.rcvTime(uIdx) / binSize) + frr.rcvOnset(uIdx);
    plot(t(idxRcv), fr(idxRcv), 'o', 'MarkerEdgeColor', 'none',...
        'MarkerFaceColor', [clrUnits(iUnit, :)], ...
        'MarkerSize', 8, 'LineWidth', 2);

    % Plot fit curve (only recovery portion)
    fitCurve = frr.fitCurve(uIdx, :);
    [~, idxRcv] = min(fitCurve);
    plot(t(idxRcv:end), fitCurve(idxRcv:end), '--', ...
        'Color', [clrUnits(iUnit, :), 1], 'LineWidth', 2);
    
    % Add asterisk at recovery time
    idxRcv = round(frr.mf.rcvTime(uIdx) / binSize) + frr.mf.rcvOnset(uIdx);
    plot(t(idxRcv), fitCurve(idxRcv), 'diamond', 'MarkerEdgeColor', 'none',...
        'MarkerFaceColor', clrUnits(iUnit, :), ...
        'MarkerSize', 8, 'LineWidth', 2);

    % Fill area of spk deficit
    tRcv = t(idxPert:end);
    frRcv = fr(idxPert:end);
    bslLine = ones(size(tRcv)) * frr.mf.frBsl(uIdx);
    tRcv = tRcv(:)';
    frRcv = frRcv(:)';
    bslLine = bslLine(:)';
    fill(hAx, [tRcv, fliplr(tRcv)], [bslLine, fliplr(frRcv)],...
        clrUnits(iUnit, :), 'FaceAlpha', 0.1, 'EdgeColor', 'none',...
        'HandleVisibility', 'off');   
end

% Mark perturbation onset
xline(pertTime, 'r--', {'Pert. Onset'}, 'LineWidth', 1,...
    'LabelOrientation', 'horizontal', 'HandleVisibility', 'off', ...
    'LabelHorizontalAlignment', 'left');

% Formatting
xlabel('Time (min)');
ylabel('Firing Rate (Hz)');
xlim([0, 24]);
xticks([0 : 6 : 24]);
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'wide', 'axHeight', 300);

