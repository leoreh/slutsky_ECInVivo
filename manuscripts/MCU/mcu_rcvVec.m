%% ========================================================================
%  MCU RECOVERY VECTOR ANALYSIS (In Vivo Pseudo-Tracking)
%  ========================================================================
%  Analyzes the recovery "strategy" of neurons in vivo by creating
%  "synthetic units" via quantile matching.
%
%  Since we cannot track the same units over days in vivo, we rank units
%  by firing rate in Baseline and BAC3 (Steady State) and pair them
%  based on their rank (Quantile Matching).
%
%  Then, we perform the same vector analysis as in MEA:
%  - Calculate a vector from Baseline (0,0) to Steady State (dSingle, dBurst)
%  - Analyze Magnitude (Strength) and Angle (Strategy)
%
%  ========================================================================



%% ========================================================================
%  LOAD DATA
%  ========================================================================

% Load table
basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];
presets = {'frNet', 'brst'};
tbl = mcu_tblVivo('basepaths', basepaths, 'flgClean', true, ...
    'presets', presets);

% Filter: RS units only
tblLme = tbl(tbl.unitType == 'RS', :);

% Filter: Only BSL and BAC3 days
idxDay = ismember(tblLme.Day, {'BSL', 'BAC3'});
tblLme = tblLme(idxDay, :);

% Assert no zero values for logs. This is instead of a pseudocount.
tblTrans = tbl_trans(tblLme, 'flg0', true, 'verbose', true);


%% ========================================================================
%  ASSERT EQUAL DISTRIBUTIONS & FILTERING
%  ========================================================================
%  Check if the data is consistent and follows expected distributions (Log-Normal).
%  Filter out low firing rate units that deviate from the distribution.

flgQq = false;
cutoff_Z = -1.8; % User defined cut-off (Visual Inspection)

if flgQq
   hFig = mcu_frQq(tbl, 'dataSet', 'vivo', 'cutoff_Z', cutoff_Z);
end

% --- FILTERING ---
cutoff_p = normcdf(cutoff_Z);
% Define threshold based on Control Baseline (canonical)
frBslRef = tblLme.fr(tblLme.Day == 'BSL' & tblLme.Group == 'Control'); 
cutoff_Hz = prctile(frBslRef, cutoff_p * 100);

fprintf('\n================================================================\n');
fprintf(' FILTERING LOW FR UNITS\n');
fprintf('================================================================\n');
fprintf('Cut-off Z-score: %.2f (%.1f%%)\n', cutoff_Z, cutoff_p*100);
fprintf('Cut-off FR     : %.4f Hz\n', cutoff_Hz);

% Apply Filter
nBefore = height(tblLme);

% Identify Units to keep (Must have BSL FR > cutoff)
goodIdx = tblLme.fr >= cutoff_Hz;

% Filter original table
tblLme = tblLme(goodIdx, :);

fprintf('Removed %d units (%.1f%%). Remaining: %d units\n', ...
    sum(~goodIdx), sum(~goodIdx)/length(goodIdx)*100, sum(goodIdx));
fprintf('================================================================\n\n');


%% ========================================================================
%  PSEUDO-TRACKING (QUANTILE MATCHING)
%  ========================================================================
%  Create 'Synthetic Units' by matching BSL and BAC3 units by rank.

nBins = 7; % Deciles

% Initialize container for synthetic table
tblSynth = match_qntl(tblLme, nBins, 'flgPool', false);

% tblGUI_scatHist(tblSynth, 'xVar', 'pBspk_BSL', 'yVar', 'fr_BAC3', 'grpVar', 'Group');


%% ========================================================================
%  VECTOR ANALYSIS
%  ========================================================================

% 1. Calculate Deltas
% -------------------
% dBrst = log(frBspk_BAC3 / frBspk_BSL)
% dSngl = log(frSspk_BAC3 / frSspk_BSL)

tblSynth.dBrst = log((tblSynth.frBspk_BAC3) ./ (tblSynth.frBspk_BSL));
tblSynth.dSngl = log((tblSynth.frSspk_BAC3) ./ (tblSynth.frSspk_BSL));

% 2. Calculate Vector Properties
% ------------------------------
[tblSynth.vecTheta, tblSynth.vecR] = cart2pol(tblSynth.dSngl, tblSynth.dBrst);
tblSynth.vecDeg = rad2deg(tblSynth.vecTheta);

% 3. Calculate Strategy Index (D)
% -------------------------------
tblSynth.D = tblSynth.dBrst - tblSynth.dSngl;

%% ========================================================================
%  PLOTTING: VECTOR STRATEGY
%  ========================================================================




