%% ========================================================================
%  ANALYSE BURSTS (IN VIVO)
%  ========================================================================

% Load table
basepaths = unique([mcu_basepaths('wt_bsl'), mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl')]);
basepaths = unique([mcu_basepaths('wt_bsl'), mcu_basepaths('mcu_bsl')]);
[tblVivo, ~, ~, ~] = mcu_tblVivo('basepaths', basepaths, 'presets', {''});

% Load states and spikes
cfgState = as_loadConfig;
vars = {'sleep_states', 'spikes'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

% Sweeping Params
isiSweep = [0.003, 0.005, 0.006, 0.008, 0.01, 0.02, 0.05, 0.1];
spkSweep = 2:5;


%% ========================================================================
%  BURST DETECTION (GRID)
%  ========================================================================
fprintf('[BRST_SWEEP] Starting Parameter Sweep (Detection)...\n');

nFiles = length(basepaths);
nIsi = length(isiSweep);
nSpk = length(spkSweep);

% Initialize Grid
burstGrid = cell(nFiles, nIsi, nSpk);

for iFile = 1:nFiles
    
    spktimes = v(iFile).spikes.times;
    
    for iIsi = 1 : nIsi
        for iSpk = 1 : nSpk
            
            minSpks = spkSweep(iSpk);
            isiStart = isiSweep(iIsi);
            isiEnd = isiStart * 2;
            minIBI = isiEnd;
            minDur = 0;
            
            % Detect
            burstGrid{iFile, iIsi, iSpk} = burst_detect(spktimes, ...
                'minSpks', minSpks, ...
                'isiStart', isiStart, ...
                'isiEnd', isiEnd, ...
                'minDur', minDur, ...
                'minIBI', minIBI, ...
                'flgForce', true, 'flgSave', false, 'flgPlot', false);
        end
    end
    fprintf('[BRST_SWEEP] Detection: File %d/%d...\n', iFile, nFiles);
end


%% ========================================================================
%  CALCULATE STATISTICS 
%  ========================================================================
fprintf('[BRST_SWEEP] Starting Parameter Sweep (Statistics)...\n');

% State
idxState = strcmp(cfgState.names, "NREM");

% -------------------------------------------------------------------------
% Pre-Process Data (Filter Spikes & Calculate Duration ONCE)
% -------------------------------------------------------------------------
fprintf('[BRST_SWEEP] Pre-processing file data (Spikes & Durations)...\n');
sData = struct('nremSpks', {}, 'nremDur', {}, 'nSpkTot', {});

for iFile = 1:nFiles
    
    spktimes = v(iFile).spikes.times;
    boutTimes = intervals(v(iFile).ss.bouts.times{idxState});
    
    % boutTimes = intervals([0 inf]);

    % Filter spikes: Keep only those within NREM bouts
    spks = cellfun(@(x) boutTimes.restrict(x, 'flgShift', false), spktimes, 'uni', false);
    
    % Calculate Duration (Sum of NREM bouts)
    nremDur = boutTimes.dur;
    
    % Store
    sData(iFile).nremSpks = spks;
    sData(iFile).nremDur  = nremDur;
    sData(iFile).nSpkTot  = cellfun(@length, spks);
end


% -------------------------------------------------------------------------
% Parameter Sweep Loop
% -------------------------------------------------------------------------
for iIsi = 1 : nIsi
    for iSpk = 1 : nSpk
        
        % Dynamic Field 
        strPrm = sprintf('s%d_i%03d', spkSweep(iSpk), round(isiSweep(iIsi)*1000));
        fPBspk = ['pBurst_', strPrm];
        
        vecPBspk = [];

        for iFile = 1:nFiles
            
            % Get Pre-processed Data
            currData  = sData(iFile);
            nremDur   = currData.nremDur;
            nSpkTot   = currData.nSpkTot; % [nUnits x 1]
            boutTimes = intervals(v(iFile).ss.bouts.times{idxState});
            
            % boutTimes = intervals([0 inf]);


            burst = burstGrid{iFile, iIsi, iSpk};
            nb   = length(burst.times);
            
            % Initialize File Stats
            pBurst  = nan(nb, 1);
            
            for iUnit = 1:nb
                % Check if bursts are fully contained in NREM
                t = burst.times{iUnit};
                if isempty(t)
                    pBurst(iUnit)  = 0;
                    continue; 
                end

                % Filter Bursts: Keep only those fully inside NREM bouts
                keepIdx = boutTimes.contains(t);

                % Count Spikes in Valid Bursts
                bSizeVals = burst.bSize{iUnit}(keepIdx);
                totBSpk   = sum(bSizeVals);
                
                % Calculate Metrics
                
                if nSpkTot(iUnit) > 0
                    pBurst(iUnit) = totBSpk / nSpkTot(iUnit);
                else
                    pBurst(iUnit) = 0;
                end
            end
             
            % Store Results (Stacking)
             % Store Results (Stacking)
            vecPBspk = [vecPBspk; pBurst];
        end
        
        % Store in Main Table
        tblVivo.(fPBspk) = vecPBspk;
    end
    fprintf('[BRST_SWEEP] Statistics: %d/%d...\n', iIsi, nIsi);
end


%% ========================================================================
%  METRICS (LME & CORRELATION)
%  ========================================================================
fprintf('[BRST_SWEEP] Starting Parameter Sweep (Metrics)...\n');

tblRes = table();

% Filter Table
idxRs = tblVivo.unitType == 'RS';
tblLme = tblVivo(idxRs, :);
idxWt = tblLme.genotype == 'Control';

for iIsi = 1 : nIsi
    for iSpk = 1 : nSpk
        
        strPrm = sprintf('s%d_i%03d', spkSweep(iSpk), round(isiSweep(iIsi)*1000));
        fPBspk = ['pBurst_', strPrm];
        
        r = struct();
        r.spkThr = spkSweep(iSpk);
        r.isiThr = isiSweep(iIsi);
        
        % LME (Group Effect)
        mdl = sprintf('%s ~ genotype + (1|sbjID)', fPBspk);
        lme = lme_analyse(tblLme, mdl, 'dist', 'logit-normal', ...
            'flgPlot', false, 'verbose', false);
        
        idxGrp = find(strncmpi(lme.Coefficients.Name, 'genotype', 8));
        if ~isempty(idxGrp)
            r.tStatGroup = lme.Coefficients.tStat(idxGrp(1));
        else
            r.tStatGroup = NaN;
        end
        r.AIC = lme.ModelCriterion.AIC;
        
        % Correlations
        r.corrFr = corr(tblLme.fr(idxWt), tblLme.(fPBspk)(idxWt), ...
            'Type', 'Spearman', 'Rows', 'complete');
        
        r.corrRoy = corr(tblLme.bRoy(idxWt), tblLme.(fPBspk)(idxWt), ...
            'Type', 'Spearman', 'Rows', 'complete');
        
        % Zero Inflation
        currPBspk = tblLme.(fPBspk)(idxWt);
        r.pZero = sum(currPBspk == 0) / length(currPBspk) * 100;

        tblRes = [tblRes; struct2table(r)];
    end
end


%% ========================================================================
%  PLOT RESULTS
%  ========================================================================

% Convert Table to Matrices for Heatmaps
matTStatGrp  = unstack(tblRes(:, {'spkThr', 'isiThr', 'tStatGroup'}), 'tStatGroup', 'spkThr');
matAIC       = unstack(tblRes(:, {'spkThr', 'isiThr', 'AIC'}), 'AIC', 'spkThr');
matCorrFr    = unstack(tblRes(:, {'spkThr', 'isiThr', 'corrFr'}), 'corrFr', 'spkThr');
matCorrRoy   = unstack(tblRes(:, {'spkThr', 'isiThr', 'corrRoy'}), 'corrRoy', 'spkThr');
mat0         = unstack(tblRes(:, {'spkThr', 'isiThr', 'pZero'}), 'pZero', 'spkThr');

% Extract matrix data (remove first col which is isiThr label)
matTStatGrp  = table2array(matTStatGrp(:, 2:end));
matAIC       = table2array(matAIC(:, 2:end));
matCorrFr    = table2array(matCorrFr(:, 2:end));
matCorrRoy   = table2array(matCorrRoy(:, 2:end));
mat0         = table2array(mat0(:, 2:end));

figure('Name', 'Burst Detection Optimization (In Vivo)', 'Color', 'w', 'Position', [100 100 1000 800]);
tiledlayout(2, 3, 'TileSpacing', 'compact');

% 1. T-Statistic (Group Effect)
nexttile;
heatmap(spkSweep, isiSweep, matTStatGrp, 'ColorMap', parula);
xlabel('Min Spikes'); ylabel('ISI Threshold (s)');
title('Difference: t-stat (Group)');

% 2. Model Fit: AIC
nexttile;
heatmap(spkSweep, isiSweep, matAIC, 'ColorMap', flipud(parula));
xlabel('Min Spikes'); ylabel('ISI Threshold (s)');
title('Model Fit: AIC (Lower is Better)');

% 3. Correlation (vs FR)
nexttile;
heatmap(spkSweep, isiSweep, matCorrFr, 'ColorMap', parula);
xlabel('Min Spikes'); ylabel('ISI Threshold (s)');
title('Correlation (pBurst vs fr)');

% 4. Correlation (vs bRoy)
nexttile;
heatmap(spkSweep, isiSweep, matCorrRoy, 'ColorMap', parula);
xlabel('Min Spikes'); ylabel('ISI Threshold (s)');
title('Correlation (pBurst vs bRoy)');

% 5. Zero Inflation
nexttile;
heatmap(spkSweep, isiSweep, mat0, 'ColorMap', parula);
xlabel('Min Spikes'); ylabel('ISI Threshold (s)');
title('Percent Zeros');





% tblGUI_bar(tblVivo, 'yVar', 'pBurst', 'xVar', 'genotype');



%% ========================================================================
%  CORRELATION: BSL BURSTINESS vs. PERTURBATION FR
%  ========================================================================
fprintf('[BRST_SWEEP] Starting Correlation Analysis (BSL pBurst vs BAC3 FR)...\n');

% -------------------------------------------------------------------------
% 1. Load Perturbation Data (wt_bac3 & mcu_bac3)
% -------------------------------------------------------------------------
bpBac3 = [mcu_basepaths('wt_bac3'), mcu_basepaths('mcu_bac3')];
[tblBac3, ~, ~, ~] = mcu_tblVivo('basepaths', bpBac3, ...
    'presets', {'frNet'}); % Load basic stats (we need FR)
    
% Filter: Only RS units
tblBac3 = tblBac3(tblBac3.unitType == 'RS', :);

% Aggregate by Mouse (Mean FR)
% Note: We want the mean FR per mouse during Baclofen Day 3
tblBac3Mouse = groupsummary(tblBac3, {'sbjID', 'genotype'}, 'mean', 'fr');
tblBac3Mouse.Properties.VariableNames{'mean_fr'} = 'frBac3';
tblBac3Mouse(:, 'GroupCount') = [];



% -------------------------------------------------------------------------
% 2. Aggregate Data (Single Table)
% -------------------------------------------------------------------------
fprintf('[BRST_SWEEP] Aggregating Baseline Data per Mouse...\n');

% Aggregate Baseline Data (All numeric variables, including pBurst columns)
% This creates a single table with one row per mouse/group
tblBslMouse = groupsummary(tblLme, {'sbjID', 'genotype'}, 'mean', vartype('numeric'));
tblBslMouse(:, 'GroupCount') = [];

% Join with Perturbation Data
% Result: Single table with 'mean_pBurst_...' and 'frBac3'
tblCorr = innerjoin(tblBslMouse, tblBac3Mouse, 'Keys', {'sbjID', 'genotype'});

% Calculate FR Recovery (% of Baseline)
tblCorr.Rcv = (tblCorr.frBac3 ./ tblCorr.mean_fr) * 100;


% -------------------------------------------------------------------------
% 3. Plotting Loop (Correlation)
% -------------------------------------------------------------------------
fprintf('[BRST_SWEEP] Plotting Correlations...\n');

figure('Name', 'Burst Sweep: BSL pBurst vs BAC3 FR Recovery', 'Color', 'w', ...
    'Position', [100 100 1200 900]);
tiledlayout(nSpk, nIsi, 'TileSpacing', 'compact', 'Padding', 'compact');

colWt = [0 0 0];
colMcu = [1 0 0];

for iSpk = 1 : nSpk
    for iIsi = 1 : nIsi
        
        % Current Parameter
        strPrm = sprintf('s%d_i%03d', spkSweep(iSpk), round(isiSweep(iIsi)*1000));
        fPBspk = ['mean_pBurst_', strPrm]; % Note: groupsummary adds 'mean_' prefix
        
        % Get Data Vectors
        pBurstData = tblCorr.(fPBspk);
        yData = tblCorr.Rcv;
        
        % INDICES (Logic is now on rows of tblCorr, i.e., mice)
        idxWt = tblCorr.genotype == 'Control';
        idxMcu = tblCorr.genotype == 'MCU-KO';
        
        % Correlations
        [rWt, pWt] = corr(pBurstData(idxWt), yData(idxWt), ...
            'Type', 'Spearman', 'Rows', 'complete');
        [rMcu, pMcu] = corr(pBurstData(idxMcu), yData(idxMcu), ...
            'Type', 'Spearman', 'Rows', 'complete');
        
        % Plot
        nexttile;
        hold on;
        
        % Scatter
        scatter(pBurstData(idxWt), yData(idxWt), 50, 'Filled', ...
            'MarkerFaceColor', colWt, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.6);
        scatter(pBurstData(idxMcu), yData(idxMcu), 50, 'Filled', ...
            'MarkerFaceColor', colMcu, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.6);
        
        % Lines
        % Note: lsline fits to all visible data in the axes, distinguishing by color property if possible
        % But standard lsline calculates regress on each line object found.
        hLines = lsline;
        if length(hLines) >= 2
            % Re-assign colors/width explicitly to match scatter
            % Handles are usually in reverse creation order
            hLines(1).Color = colMcu; hLines(1).LineWidth = 1.5; 
            hLines(2).Color = colWt;  hLines(2).LineWidth = 1.5; 
        end
        
        axis square;
        grid on;
        
        % Formatting
        title(sprintf('S:%d, I:%.3f\nWT: r=%.2f | MCU: r=%.2f', ...
            spkSweep(iSpk), isiSweep(iIsi), rWt, rMcu), ...
            'FontWeight', 'normal', 'FontSize', 9);
        
        if iSpk == nSpk
            xlabel('BSL pBurst');
        end
        if iIsi == 1
            ylabel('FR Recovery (%)');
        end
        
        set(gca, 'TickDir', 'out', 'Box', 'off');
        hold off;
    end
end

fprintf('[BRST_SWEEP] Correlation analysis complete.\n');


% tblGUI_scatHist(tblCorr, 'xVar', 'pBurst', 'yVar', 'funcon_fish', 'grpVar', 'genotype');
