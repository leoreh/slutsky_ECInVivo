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
dataSet = 'mea';
dataSet = 'vivo';

if flgLoad
    if strcmp(dataSet, 'mea')
        presets = {'steadyState'};
        [tblMea, xVec, basepaths, v] = mcu_tblMea('presets', presets, 'flgOtl', true);

    elseif strcmp(dataSet, 'vivo')
        basepaths = [mcu_basepaths('wt_bsl'), mcu_basepaths('mcu_bsl'), ...
            mcu_basepaths('wt_bac3'), mcu_basepaths('mcu_bac3')];
        basepaths(contains(basepaths, 'lh137')) = [];
        presets = {'burst'};
        tblVivo = mcu_tblVivo('basepaths', basepaths, 'flgClean', true, ...
            'presets', presets);
        tblVivo.day(tblVivo.day == "BAC_ON") = "BAC3";
        tblVivo.day = removecats(tblVivo.day, {'BAC_ON'});
        
        % Filter RS
        tblVivo = tblVivo(tblVivo.unitType == "RS", :);
        tblVivo.unitType = [];

        % Assert no zero values (instead of a pseudocount)
        tblVivo = tbl_trans(tblVivo, 'flg0', true, 'verbose', true);
        
    end
end

if strcmp(dataSet, 'mea')
    tblRaw = tblMea;

elseif strcmp(dataSet, 'vivo')
    tblRaw = tblVivo;
end

%% ========================================================================
%  ASSERT EQUAL DISTRIBUTIONS
%  ========================================================================

flgPool = false;

if strcmp(dataSet, 'mea')
    cutoff_z = -2;
    nBins = 5;
elseif strcmp(dataSet, 'vivo')
    cutoff_z = -1.0;
    nBins = 4;
end

[hFig, qqData] = mcu_rcvQq(tblRaw, 'var', 'fr', 'dataSet', dataSet, 'cutoff_z', cutoff_z);


% Define threshold based on Control Baseline
frRef = log(tblRaw.fr(tblRaw.genotype == 'Control'));
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
    varsTbl = {'fr', 'frBurst', 'frSingle', 'pBurst'};
    tBsl = tblRaw(:, [{'genotype', 'sbjID'}, varsTbl]);
    tBsl.day = repmat({'BSL'}, height(tBsl), 1);
    tBsl.day = categorical(tBsl.day);

    % Steady State Table (Rename ss_ vars)
    varsSs = {'ss_fr', 'ss_frBurst', 'ss_frSingle', 'ss_pBurst'};
    tSs = tblRaw(:, [{'genotype', 'sbjID'}, varsSs]);
    tSs = renamevars(tSs, {'ss_fr', 'ss_frBurst', 'ss_frSingle', 'ss_pBurst'}, ...
        {'fr',    'frBurst',    'frSingle',    'pBurst'});
    tSs.day = repmat({'BAC3'}, height(tSs), 1); % match_qntl hardcodes 'BAC3' as the 2nd day
    tSs.day = categorical(tSs.day);

    % Combine
    tblRaw = [tBsl; tSs];
end

%% ========================================================================
%  PSEUDO-TRACKING
%  ========================================================================
fprintf('VALIDATION: Running match_qntl (nBins=%d)...\n', nBins);

tbl = mcu_matchQntl(tblRaw, nBins, 'flgPool', flgPool, ...
    'var', 'fr', 'avgType', 'geomean');


%% ========================================================================
%  ADD CALCS
%  ========================================================================

% logit pBurst
tblTrans = tbl_trans(tbl, 'varsInc', {'pBurst', 'ss_pBurst'}, 'logBase', 'logit');
tbl.pBurst_trans = tblTrans.pBurst;
tbl.ss_pBurst_trans = tblTrans.ss_pBurst;

% Relative
tbl.dBrst_rel = log((tbl.ss_frBurst) ./ (tbl.frBurst));
tbl.dSngl_rel = log((tbl.ss_frSingle) ./ (tbl.frSingle));
tbl.dFr = log((tbl.ss_fr) ./ (tbl.fr));
tbl.dpBurst = (tbl.ss_pBurst_trans) - (tbl.pBurst_trans);

% Absolute
tbl.dBrst_abs = (tbl.ss_frBurst - tbl.frBurst);
tbl.dSngl_abs = (tbl.ss_frSingle - tbl.frSingle);
tbl.dFr_abs = (tbl.ss_fr - tbl.fr);


%% ========================================================================
%  PLAY
%  ========================================================================
% Repeat regression analyses

hFig = mcu_rcvSpace(tbl);

% tblGUI_scatHist(tbl, 'xVar', 'pBurst', 'yVar', 'ss_frBurst', 'grpVar', 'genotype');
