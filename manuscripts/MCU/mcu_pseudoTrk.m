%% ========================================================================
%  MCU PSEUDO-TRACKING
%  ========================================================================
% Creates pairs of synthetic "neurons" via quantile matching to validate
% the regression analyses.
%
% Since we cannot track the same units over days in vivo, we rank units
% by firing rate in Baseline and BAC3 (Steady State) and pair them
% based on their rank (Quantile Matching).
%
% Designed for in vivo but can work on the mea data set as well.

%% ========================================================================
%  LOAD DATA
%  ========================================================================

flgLoad = false;
dataSet = 'vivo';

if flgLoad
    if strcmp(dataSet, 'mea')
        presets = {'steadyState'};
        [tblMea, xVec, basepaths, v] = mcu_tblMea('presets', presets, 'flgOtl', true);

    elseif strcmp(dataSet, 'vivo')
        basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];
        presets = {'frNet', 'brst'};
        tblVivo = mcu_tblVivo('basepaths', basepaths, 'flgClean', true, ...
            'presets', presets);  
    end
end

% Prep table
if strcmp(dataSet, 'mea')
    tblRaw = tblMea;

elseif strcmp(dataSet, 'vivo')
    % Filter
    tblRaw = tblVivo(tblVivo.unitType == 'RS', :);

    idxDay = ismember(tblRaw.Day, {'BSL', 'BAC3'});
    tblRaw = tblRaw(idxDay, :);

    % Assert no zero values for logs (instead of a pseudocount)
    tblRaw = tbl_trans(tblRaw, 'flg0', true, 'verbose', true);
end

%% ========================================================================
%  ASSERT EQUAL DISTRIBUTIONS
%  ========================================================================

flgPool = false;

if strcmp(dataSet, 'mea')
    cutoff_z = -2;
    nBins = 10;
elseif strcmp(dataSet, 'vivo')
    cutoff_z = -1;
    nBins = 3;
end

hFig = mcu_rcvQq(tblRaw, 'var', 'frSspk', 'dataSet', dataSet, 'cutoff_z', cutoff_z);

% Define threshold based on Control Baseline
frRef = log(tblRaw.fr(tblRaw.Group == 'Control'));
mu = mean(frRef);
sigma = std(frRef);
cutoff_Hz = exp(mu + cutoff_z * sigma);
cutoff_p = normcdf(cutoff_z); 

% Filter
nBefore = height(tblRaw);
goodIdx = tblRaw.fr >= cutoff_Hz;
tblRaw = tblRaw(goodIdx, :);

fprintf('\n================================================================\n');
fprintf(' FILTERING LOW FR UNITS\n');
fprintf('================================================================\n');
fprintf('Cut-off Z-score: %.2f (%.1f%%)\n', cutoff_z, cutoff_p*100);
fprintf('Cut-off FR     : %.4f Hz\n', cutoff_Hz);
fprintf('Removed %d units (%.1f%%). Remaining: %d units\n', ...
    sum(~goodIdx), sum(~goodIdx)/length(goodIdx)*100, sum(goodIdx));
fprintf('================================================================\n\n');


% --- CONVERT TO LONG FORMAT ---
if strcmp(dataSet, 'mea')
    % Baseline Table
    varsTbl = {'fr', 'frBspk', 'frSspk', 'pBspk', 'pBspk_trans'};
    tBsl = tblRaw(:, [{'Group', 'Name'}, varsTbl]);
    tBsl.Day = repmat({'BSL'}, height(tBsl), 1);
    tBsl.Day = categorical(tBsl.Day);

    % Steady State Table (Rename ss_ vars)
    varsSs = {'frSs', 'ss_frBspk', 'ss_frSspk', 'ss_pBspk', 'ss_pBspk_trans'};
    tSs = tblRaw(:, [{'Group', 'Name'}, varsSs]);
    tSs = renamevars(tSs, {'frSs', 'ss_frBspk', 'ss_frSspk', 'ss_pBspk', 'ss_pBspk_trans'}, ...
        {'fr',    'frBspk',    'frSspk',    'pBspk', 'pBspk_trans'});
    tSs.Day = repmat({'BAC3'}, height(tSs), 1); % match_qntl hardcodes 'BAC3' as the 2nd day
    tSs.Day = categorical(tSs.Day);

    % Combine
    tblRaw = [tBsl; tSs];
end

%% ========================================================================
%  PSEUDO-TRACKING
%  ========================================================================
fprintf('VALIDATION: Running match_qntl (nBins=%d)...\n', nBins);

tbl = mcu_matchQntl(tblRaw, nBins, 'flgPool', flgPool, ...
    'var', 'fr', 'avgType', 'mean');


%% ========================================================================
%  ADD CALCS
%  ========================================================================

% Naming convention
tbl.frSs = tbl.ss_fr;
removevars(tbl, 'ss_fr');

% logit pBspk
tblTrans = tbl_trans(tbl, 'varsInc', {'pBspk', 'ss_pBspk'}, 'logBase', 'logit');
tbl.pBspk_trans = tblTrans.pBspk;
tbl.ss_pBspk_trans = tblTrans.ss_pBspk;

% Relative
tbl.dBrst_rel = log((tbl.ss_frBspk) ./ (tbl.frBspk));
tbl.dSngl_rel = log((tbl.ss_frSspk) ./ (tbl.frSspk));
tbl.dFr = log((tbl.frSs) ./ (tbl.fr));
tbl.dpBspk = (tbl.ss_pBspk_trans) - (tbl.pBspk_trans);

% Absolute
tbl.dBrst_abs = (tbl.ss_frBspk - tbl.frBspk);
tbl.dSngl_abs = (tbl.ss_frSspk - tbl.frSspk);
tbl.dFr_abs = (tbl.frSs - tbl.fr);


%% ========================================================================
%  PLAY
%  ========================================================================
% Repeat regression analyses

hFig = mcu_rcvSpace(tbl);
