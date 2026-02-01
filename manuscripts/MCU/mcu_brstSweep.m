%% ========================================================================
%  ANALYSE BURSTS (IN VIVO)
%  ========================================================================

% Load table
basepaths = unique([mcu_basepaths('wt_bsl'), mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl')]);
[tblVivo, ~, ~, ~] = mcu_tblVivo('basepaths', basepaths, 'presets', {''});

% Load states and spikes
cfgState = as_loadConfig;
vars = {'sleep_states', 'spikes'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

% State
idxState = contains(cfgState.names, "NREM");

% Sweeping Params
isiSweep = [0.005, 0.01, 0.02, 0.05, 0.1];
spkSweep = 2:5;


%% ========================================================================
%  BURST DETECTION (GRID)
%  ========================================================================
fprintf('[BRST_SWEEP] Starting Parameter Sweep (Detection)...\n');

nFiles = length(basepaths);
nIsi = length(isiSweep);
nSpk = length(spkSweep);

% Initialize Grid
brstGrid = cell(nFiles, nIsi, nSpk);

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
            brstGrid{iFile, iIsi, iSpk} = brst_detect(spktimes, ...
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

% -------------------------------------------------------------------------
% Pre-Process Data (Filter Spikes & Calculate Duration ONCE)
% -------------------------------------------------------------------------
fprintf('[BRST_SWEEP] Pre-processing file data (Spikes & Durations)...\n');
sData = struct('nremSpks', {}, 'nremDur', {}, 'nSpkTot', {});

for iFile = 1:nFiles
    
    spktimes = v(iFile).spikes.times;
    boutTimes = ints(v(iFile).ss.bouts.times{idxState});
    
    % Filter spikes: Keep only those within NREM bouts
    tic
    spks = cellfun(@(x) RestrictInts(x, boutTimes.list), spktimes, 'uni', false);
    toc

    tic
    spks2 = cellfun(@(x) boutTimes.restrict(x), spktimes, 'uni', false);
    toc
    

    
    % Calculate Duration (Sum of NREM bouts)
    if isempty(boutTimes)
        nremDur = 0;
    else
        nremDur = sum(boutTimes(:,2) - boutTimes(:,1));
    end
    
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
        fSib = ['sib_', strPrm];
        fFrB = ['frB_', strPrm];
        fFrS = ['frS_', strPrm];
        
        vecSib = [];
        vecFrB = [];
        vecFrS = [];

        for iFile = 1:nFiles
            
            % Get Pre-processed Data
            currData  = sData(iFile);
            nremDur   = currData.nremDur;
            nSpkTot   = currData.nSpkTot; % [nUnits x 1]
            boutTimes = v(iFile).ss.bouts.times{idxState};
            
            brst = brstGrid{iFile, iIsi, iSpk};
            nb   = length(brst.times);
            
            % Initialize File Stats
            pBspk  = nan(nb, 1);
            frBspk = zeros(nb, 1);
            frSspk = zeros(nb, 1);
            
            for iUnit = 1:nb
                % Check if bursts are fully contained in NREM
                t = brst.times{iUnit};
                if isempty(t)
                    pBspk(iUnit)  = 0;
                    frBspk(iUnit) = 0;
                    if nremDur > 0
                        frSspk(iUnit) = nSpkTot(iUnit) / nremDur;
                    end
                    continue; 
                end

                % Filter Bursts: Keep only those fully inside NREM bouts
                keepIdx = InIntervals(t(:,1), boutTimes) & InIntervals(t(:,2), boutTimes);
                
                % Count Spikes in Valid Bursts
                nBspkVals = brst.nBspk{iUnit}(keepIdx);
                totBSpk   = sum(nBspkVals);
                
                % Calculate Metrics
                if nremDur > 0
                    frBspk(iUnit) = totBSpk / nremDur;
                    frSspk(iUnit) = (nSpkTot(iUnit) - totBSpk) / nremDur;
                end
                
                if nSpkTot(iUnit) > 0
                    pBspk(iUnit) = totBSpk / nSpkTot(iUnit);
                else
                    pBspk(iUnit) = 0;
                end
            end
             
            % Store Results (Stacking)
            vecSib = [vecSib; pBspk];
            vecFrB = [vecFrB; frBspk];
            vecFrS = [vecFrS; frSspk];
        end
        
        % Store in Main Table
        tblVivo.(fSib) = vecSib;
        tblVivo.(fFrB) = vecFrB;
        tblVivo.(fFrS) = vecFrS;
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
idxWt = tblLme.Group == 'Control';

for iIsi = 1 : nIsi
    for iSpk = 1 : nSpk
        
        strPrm = sprintf('s%d_i%03d', spkSweep(iSpk), round(isiSweep(iIsi)*1000));
        fSib = ['sib_', strPrm];
        
        r = struct();
        r.spkThr = spkSweep(iSpk);
        r.isiThr = isiSweep(iIsi);
        
        % LME (Group Effect)
        mdl = sprintf('%s ~ Group + (1|Name)', fSib);
        lme = lme_analyse(tblLme, mdl, 'dist', 'logit-normal', ...
            'flgPlot', false, 'verbose', false);
        
        idxGrp = find(strncmpi(lme.Coefficients.Name, 'Group', 5));
        if ~isempty(idxGrp)
            r.tStatGroup = lme.Coefficients.tStat(idxGrp(1));
        else
            r.tStatGroup = NaN;
        end
        r.AIC = lme.ModelCriterion.AIC;
        
        % Correlations
        r.corrFr = corr(tblLme.fr(idxWt), tblLme.(fSib)(idxWt), ...
            'Type', 'Spearman', 'Rows', 'complete');
        
        if ismember('bRoy', tblLme.Properties.VariableNames)
            r.corrRoy = corr(tblLme.bRoy(idxWt), tblLme.(fSib)(idxWt), ...
                'Type', 'Spearman', 'Rows', 'complete');
        else
            r.corrRoy = NaN;
        end
        
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

% Extract matrix data (remove first col which is isiThr label)
matTStatGrp  = table2array(matTStatGrp(:, 2:end));
matAIC       = table2array(matAIC(:, 2:end));
matCorrFr    = table2array(matCorrFr(:, 2:end));
matCorrRoy   = table2array(matCorrRoy(:, 2:end));

figure('Name', 'Burst Detection Optimization (In Vivo)', 'Color', 'w', 'Position', [100 100 1000 800]);
tiledlayout(2, 2, 'TileSpacing', 'compact');

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
title('Correlation (sib vs fr)');

% 4. Correlation (vs bRoy)
nexttile;
heatmap(spkSweep, isiSweep, matCorrRoy, 'ColorMap', parula);
xlabel('Min Spikes'); ylabel('ISI Threshold (s)');
title('Correlation (sib vs bRoy)');




tblGUI_bar(tblLme, 'yVar', 'pBspk', 'xVar', 'Group');
