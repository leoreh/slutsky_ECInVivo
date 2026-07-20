



%% ========================================================================
%  CAG-MCU-KO - METRICS (st & bursts)
%  ========================================================================
% Five MCU-KO recordings (CAG cohort, 'ra' path key; TDT, fs ~24.4 kHz).
% Spike timing metrics, bursts and burst stats, computed with the same params
% as the wt / mcu pipeline (see mcu_wrapper). Unit classification is handled
% separately in the next section.

basepaths = mcu_basepaths('ra');
nFiles = length(basepaths);

% vars
vars = {'spikes'};

% load state vars
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

% Burst detection params (identical to mcu_wrapper)
isiStart = 0.006;
minSpks  = 2;
isiEnd   = isiStart * 2;
minIBI   = isiEnd;
minDur   = 0;

for iFile = 1 : nFiles

    % file
    basepath = basepaths{iFile};
    [~, basename] = fileparts(basepath);
    cd(basepath)

    session = CE_sessionTemplate(pwd, 'viaGUI', false,...
        'forceDef', true, 'forceL', true, 'saveVar', true);
    basepath = session.general.basePath;
    nchans = session.extracellular.nChannels;
    fs = session.extracellular.sr;
    spkgrp = session.extracellular.spikeGroups.channels;
    [~, basename] = fileparts(basepath);
    
    % % Sleep Signal
    % sSig = as_prepSig([basename, '.lfp'], [basename, '.emg.dat'],...
    %     'eegCh', [1 : 4], 'emgCh', 1, 'saveVar', true, 'emgNchans', 4, 'eegNchans', 12,...
    %     'inspectSig', false, 'forceLoad', true, 'eegFs', 1250, 'emgFs', 3051.76,...
    %     'emgCf', [80 450]);
    % 
    % labelsmanfile = [basename, '.sleep_labelsMan.mat'];
    % AccuSleep_viewer(sSig, labels, labelsmanfile)
    % 
    % % classify with a network
    % netfile = 'D:\Code\slutsky_ECInVivo\lfp\SleepStates\AccuSleep\trainedNetworks\net_230212_103132.mat';
    % calData = [];
    % ss = as_classify(sSig, 'basepath', basepath, 'inspectLabels', false,...
    %     'saveVar', true, 'forceA', true, 'netfile', netfile,...
    %     'graphics', true, 'calData', calData);
    % 
    % % spktimes
    % spktimes = v(iFile).spikes.times;
    % 
    % % Spike timing metrics
    % st = spktimes_metrics('spktimes', spktimes, 'sunits', [], ...
    %     'bins', {[0, Inf]}, 'flgForce', true, 'flgSave', true, ...
    %     'flgAll', false);
    % 
    % % Burst detection
    % burst = burst_detect(spktimes, ...
    %     'minSpks', minSpks, ...
    %     'isiStart', isiStart, ...
    %     'isiEnd', isiEnd * 2, ...
    %     'minDur', minDur, ...
    %     'minIBI', minIBI, ...
    %     'flgForce', true, 'flgSave', true, 'flgPlot', false);
    % 
    % % Burst statistics
    % stats = burst_stats(burst, spktimes, 'winCalc', [], 'flgSave', true);
    
    % spike wave metrics
    swv = spkwv_metrics('basepath', basepath, 'flgSave', true,...
        'flgForce', true);

    % % Ripples
    % load([basename, '.sleep_sig.mat'], 'info');
    % ripp = ripp_wrapper('basepath', pwd, ...
    %     'win', [0 12] * 3600, ...
    %     'rippCh', info.eegCh, ...
    %     'flgPlot', true, ...
    %     'flgSave', true, ...
    %     'flgNS', true, ...
    %     'flgForce', true);
    % 
    % % Epileptiform discharges
    % ed = ed_wrapper('basepath', basepath, ...
    %     'flgSave', true, ...
    %     'flgPlot', false, ...
    %     'flgForce', true);  


end


guiPath(pwd, 'preset', 'ripp');



% Unit Class: GMM classification, then manual curation in the GUI
fetSelect = {'Asym', 'Hpk', 'TP'};
tblUnit = utypes_classify('basepaths', basepaths, 'fetSelect', fetSelect, ...
    'rsPrior', 0.97, 'regVal', 0.01, 'flgPlot', true);



%% ========================================================================
%  CAG-MCU-KO - BURSTINESS COMPARISON
%  ========================================================================
% Compare baseline burstiness of the CAG-MCU-KO mice against the wt and mcu
% baseline groups. The cohort is registered in mcu_cfg (cfg.miceCAG), so
% mcu_tblVivo labels the genotypes itself - no post hoc addcats.

basepaths = mcu_basepaths('bsl3');
tblVivo = mcu_tblVivo('basepaths', basepaths, 'presets', {'burst'}, ...
    'flgClean', true);
tblTrans = tbl_trans(tblLme, 'varsInc', {'pBurst'}, 'logBase', 'logit');

% GUI
guiTbl_bar(tblVivo, 'yVar', 'pBurst', 'xVar', 'genotype');

% LME
frml = 'br ~ genotype * fr + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblVivo, frml, 'dist', 'gamma');

frml = 'bSize ~ genotype * fr + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblVivo, frml, 'dist', 'log-normal');

frml = 'pBurst ~ genotype * fr + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblVivo, frml, 'dist', 'logit-normal');

frml = 'fr ~ genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblVivo, frml, 'dist', 'gamma');

cn = lmeMdl.CoefficientNames;
L  = zeros(1, numel(cn));
L(strcmp(cn, 'genotype_CAG-MCU-KO')) =  1;
L(strcmp(cn, 'genotype_MCU-KO'))     = -1;
est = L * lmeMdl.Coefficients.Estimate;
se  = sqrt(L * lmeMdl.CoefficientCovariance * L');
[pVal, ~, ~, df] = coefTest(lmeMdl, L, 0);

% 2 Prism
tblVivo.fr(tblVivo.genotype == 'CAG-MCU-KO')

% --- States

tblVivo = mcu_tblVivo('basepaths', basepaths, 'presets', {'spkStates'}, ...
    'flgClean', true);

% GUI
guiTbl_bar(tblVivo, 'yVar', 'pBurst', 'xVar', 'genotype', 'grpVar', 'state');


% --- 2 Prism

% One metric across states x genotypes, summarized per cell. Laid out for a
% Prism Grouped table with the "Mean, SD, N" entry format: rows are the
% states (x-axis), each genotype takes three consecutive subcolumns. Only
% the numeric block is copied - pasting the labels with it would shift
% every value one row and one column.
varPrism = 'fr';
sstates = {'WAKE', 'NREM', 'REM'};
uGeno = categories(removecats(tblVivo.genotype));
nGeno = length(uGeno);

% REM bouts are short, so a unit can enter a state row on a handful of
% spikes. Raise to gate those out
minSpks = 0;

prismMat = nan(length(sstates), nGeno * 3);
for iState = 1 : length(sstates)
    for iGeno = 1 : nGeno

        idx = tblVivo.state == sstates{iState} & ...
            tblVivo.genotype == uGeno{iGeno} & ...
            tblVivo.nSpks >= minSpks;
        dat = tblVivo.(varPrism)(idx);
        dat = dat(~isnan(dat));

        prismMat(iState, (iGeno - 1) * 3 + (1 : 3)) = ...
            [mean(dat), std(dat), numel(dat)];
    end
end

% Inspect column order
colLbls = strcat(repelem(string(uGeno(:))', 1, 3), ...
    repmat(["_mean", "_sd", "_n"], 1, nGeno));
prismTbl = array2table(prismMat, 'RowNames', sstates, ...
    'VariableNames', matlab.lang.makeValidName(colLbls));
disp(prismTbl)

% Clipboard
clipboard('copy', sprintf([repmat('%.6g\t', 1, nGeno * 3 - 1), '%.6g\n'], ...
    prismMat'));



%% ========================================================================
%  CAG-MCU-KO - SPIKE TIMING TRACES (ACG & ISI)
%  ========================================================================
% Population average of two per-unit traces, across RS units, for the three
% genotype baseline set: the narrow autocorrelogram (from st_metrics) and
% the log-ISI distribution (computed here - spktimes_metrics keeps isi
% scalars but not the histogram). Each unit contributes one normalised
% trace, so a fast unit does not dominate the mean.
%
% Both export to a Prism XY table with the "Mean, SD, N" entry format:
% column 1 is X, then three subcolumns per genotype.

basepaths = mcu_basepaths('bsl3');
[tblTrc, ~, ~, xVec] = mcu_tblVivo('basepaths', basepaths, ...
    'presets', {'acg', 'spktimes'}, 'flgClean', true);

% acg lag axis [ms]. spktimes_metrics returns nan for a unit under its 100
% spike floor and for lags no spike could reach, so those drop out per bin
xAcg = xVec.narrow(:)';


% --- ISI histogram

% Log spaced bins, 1 ms to 100 s. Each unit is normalised by its own
% in-range isi count and then by the bin width, so the trace is a density
% over log10(isi) - it integrates to 1 per decade and its height does not
% move if dLog changes. Y units are probability per decade; a plain
% fraction per bin would only be readable next to the bin width. Whole
% recording, so there are no bout edges to respect (cf. spktimes_metrics,
% which is segment aware)
dLog = 0.05;                                    % decade per bin
isiEdges = 10 .^ (-3 : dLog : 2);
xIsi = sqrt(isiEdges(1 : end - 1) .* isiEdges(2 : end)) * 1000;  % [ms]

% same exposure floor spktimes_metrics applies to the acg
minSpks = 100;

isiHist = nan(height(tblTrc), length(xIsi));
for iUnit = 1 : height(tblTrc)

    isi = diff(sort(tblTrc.spktimes{iUnit}(:)));
    if length(isi) < minSpks
        continue
    end

    cnt = histcounts(isi, isiEdges);
    isiHist(iUnit, :) = cnt / sum(cnt) / dLog;
end
tblTrc.isiHist = isiHist;
% 
% % GUI. one viewer each - guiTbl_xy matches y vars by length against xVec
% guiTbl_xy(xAcg, tblTrc, 'yVar', 'acgNarrow', 'grpVar', 'genotype', ...
%     'xLbl', 'Lag [ms]');
% guiTbl_xy(xIsi, tblTrc, 'yVar', 'isiHist', 'grpVar', 'genotype', ...
%     'xLbl', 'ISI [ms]');
% 

% --- 2 Prism

% Rows are x bins, columns are x then genotype x (mean, sd, n). n is counted
% per bin so it reports the units that actually contributed there, which the
% single GroupCount of a groupsummary would overstate
uGeno = categories(removecats(tblTrc.genotype));
nGeno = length(uGeno);

trcVars = {'acgNarrow', 'isiHist'};
trcX = {xAcg, xIsi};

prism = struct();
for iTrc = 1 : length(trcVars)

    dat = tblTrc.(trcVars{iTrc});
    mat = nan(size(dat, 2), nGeno * 3);

    for iGeno = 1 : nGeno
        datGeno = dat(tblTrc.genotype == uGeno{iGeno}, :);
        mat(:, (iGeno - 1) * 3 + (1 : 3)) = ...
            [mean(datGeno, 1, 'omitnan')', std(datGeno, 0, 1, 'omitnan')', ...
            sum(~isnan(datGeno), 1)'];
    end

    prism.(trcVars{iTrc}) = [trcX{iTrc}(:), mat];
end

% Column order, for the paste
disp(['X, ', strjoin(strcat(repelem(string(uGeno(:))', 1, 3), ...
    repmat([" mean", " sd", " n"], 1, nGeno)), ', ')])

% Clipboard. swap the field for the other trace
prismMat = prism.isiHist;
clipboard('copy', sprintf([repmat('%.6g\t', 1, size(prismMat, 2) - 1), ...
    '%.6g\n'], prismMat'));