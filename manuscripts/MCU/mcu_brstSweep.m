
%% ========================================================================
%  MEA ANALYSE BURSTS
%  ========================================================================

% Load table
basepaths = [mcu_basepaths('mea_bac')];
presets = {'spktimes', 'rcv', 'frNet'};
[tblFull, ~, ~, v] = mcu_tblMea('basepaths', basepaths, 'presets', presets([1, 3]));

% Prepare subtable of what's needed
tbl = tblFull(:, {'Name', 'Group', 'UnitID', 'fr', 'frSs', 'frTrough', ...
    'uPert', 'spktimes', 'funcon'}); 

% Clean units that were not perturbed
tbl(~tbl.uPert, :) = [];
tbl = removevars(tbl, 'uPert');

% Spike times from all sessions
spktimes = tbl.spktimes;

% Fixed params
binSize = 60;

% Baseline Window
rcv = catfields([v(:).rcv], 1);
winBsl = rcv.info.winBsl;
winBsl = [min(winBsl(:)), min(winBsl(:, 2))];

% Sweeping Params
isiVal = 0.05;

for iSweep = 1 : length(isiVal)
    
    minSpks = 3;
    isiThr = isiVal(iSweep);

    % Burst detection
    brst = brst_detect(spktimes, ...
        'minSpks', minSpks, ...
        'isiStart', isiThr, ...
        'isiEnd', isiThr * 2, ...
        'minDur', 0.005, ...
        'minIBI', 0.1, ...
        'flgForce', true, 'flgSave', false, 'flgPlot', false);

    % Burst statistics
    stats = brst_stats(brst, spktimes, 'winCalc', winBsl, ...
        'flgSave', false);
    
    tbl.sib = stats.pBspk;
end


%% ========================================================================
%  MEA PREDICTIVE POWER
%  ========================================================================

frml = 'frSs ~ sib + (1|Name)';

for iSweep = 1 : length(isiVal)

    frml = ['frSs ~ ', 'sib', ' + (1|Name)'];

    lmeMdl = lme_analyse(tbl, frml, 'dist', 'log-normal');
    fxdEffects = lmeMdl.Coefficients;
    tStat(iSweep) = fxdEffects.tStat(2);

end


%% ========================================================================
%  Full Model
%  ========================================================================

frml = 'frSs ~ (fr + sib + frTrough + funcon) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, 'dist', 'log-normal');








%% ========================================================================
%  MEA ISI VALLEY
%  ========================================================================
% Completely useless. No bi-modality. Also, optimization via predictive
% power is preferred.

presets = {'spktimes'};
[tbl, xVec, basepaths, v] = mcu_tblMea('presets', presets([1]));

idxWt = tbl.Group == 'Control';
spktimes = tbl.spktimes(idxWt);

% ISI VALLEY AS THRESHOLD
isiVal = brst_isiValley(spktimes, 'nSpks', 3);




%% ========================================================================
%  IN VIVO
%  ========================================================================

basepaths = [mcu_basepaths('wt_bsl')];
[tbl, ~, ~, ~] = mcu_tblVivo('basepaths', basepaths, 'presets', {'spktimes'});

idxUnit = tbl.UnitType == 'RS';
spktimes = tbl.spktimes(idxUnit);

% ISI VALLEY AS THRESHOLD
isiVal = brst_isiValley(spktimes, 'nSpks', 2);