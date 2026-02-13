%% ========================================================================
%  MEA RECOVERY VECTOR ANALYSIS (State Space Trajectory)
%  ========================================================================
%  Analyzes the recovery "strategy" of neurons by decomposing firing rate
%  recovery into Single-Spike and Burst-Spike components.

%#ok<*NASGU>

%% ========================================================================
%  LOAD DATA
%  ========================================================================

% Load with steady state variables
presets = {'steadyState'};
[tbl, xVec, basepaths, v] = mcu_tblMea('presets', presets, 'flgOtl', true);

% Add logit pBspk
tblTrans = tbl_trans(tbl, 'varsInc', {'pBspk', 'ss_pBspk'}, 'logBase', 'logit');
tbl.pBspk_trans = tblTrans.pBspk;
tbl.ss_pBspk_trans = tblTrans.ss_pBspk;


%% ========================================================================
%  PSEUDO-TRACKING
%  ========================================================================
% Validate "Pseudo-Tracking" (Quantile Matching) on MEA data

flgPseudo = true; 
if flgPseudo
    
    nBins = 10;

    % --- ASSERT EQUAL DISTRIBUTIONS ---
    cutoff_z = -2; % User defined cut-off (Visual Inspection)
    % cutoff_z = -6;

    hFig = mcu_rcvQq(tbl, 'var', 'frSspk', 'dataSet', 'mea', 'cutoff_z', cutoff_z);

    % Define threshold based on Control Baseline
    frRef = log(tbl.fr(tbl.Group == 'Control'));
    mu = mean(frRef);
    sigma = std(frRef);
    cutoff_Hz = exp(mu + cutoff_z * sigma);

    cutoff_p = normcdf(cutoff_z); % For display only

    % Apply Filter
    nBefore = height(tbl);

    % Identify Units to keep (Must have BSL FR > cutoff)
    goodIdx = tbl.fr >= cutoff_Hz;

    % Filter original table
    tbl = tbl(goodIdx, :);

    fprintf('\n================================================================\n');
    fprintf(' FILTERING LOW FR UNITS\n');
    fprintf('================================================================\n');
    fprintf('Cut-off Z-score: %.2f (%.1f%%)\n', cutoff_z, cutoff_p*100);
    fprintf('Cut-off FR     : %.4f Hz\n', cutoff_Hz);
    fprintf('Removed %d units (%.1f%%). Remaining: %d units\n', ...
        sum(~goodIdx), sum(~goodIdx)/length(goodIdx)*100, sum(goodIdx));
    fprintf('================================================================\n\n');


    % --- CONVERT TO LONG FORMAT ---
    
    % Baseline Table
    varsTbl = {'fr', 'frBspk', 'frSspk', 'pBspk', 'pBspk_trans'};
    tBsl = tbl(:, [{'Group', 'Name'}, varsTbl]);
    tBsl.Day = repmat({'BSL'}, height(tBsl), 1);
    tBsl.Day = categorical(tBsl.Day);

    % Steady State Table (Rename ss_ vars)
    varsSs = {'frSs', 'ss_frBspk', 'ss_frSspk', 'ss_pBspk', 'ss_pBspk_trans'}; 
    tSs = tbl(:, [{'Group', 'Name'}, varsSs]);
    tSs = renamevars(tSs, {'frSs', 'ss_frBspk', 'ss_frSspk', 'ss_pBspk', 'ss_pBspk_trans'}, ...
                           {'fr',    'frBspk',    'frSspk',    'pBspk', 'pBspk_trans'});
    tSs.Day = repmat({'BAC3'}, height(tSs), 1); % match_qntl hardcodes 'BAC3' as the 2nd day
    tSs.Day = categorical(tSs.Day);

    % Combine
    tblLong = [tBsl; tSs];
    
    % --- RUN MATCHING ---
    fprintf('VALIDATION: Running match_qntl on MEA data (nBins=%d)...\n', nBins);
    tbl = match_qntl(tblLong, nBins, 'flgPool', false, ...
        'var', 'fr', 'avgType', 'mean');
    
    tbl.frSs = tbl.ss_fr;
    fprintf('VALIDATION: Replaced real units with %d Synthetic Units.\n', height(tbl));
end

%% ========================================================================
%  PLAY 
%  ========================================================================
% Repeat regression analyses