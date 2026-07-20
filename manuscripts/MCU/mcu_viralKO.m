



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

frml = 'pBurst ~ genotype * fr + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblVivo, frml, 'dist', 'logit-normal');

% Compare (switch yVar in the GUI for frBurst, br, bSize, etc.)
guiTbl_bar(tblVivo, 'yVar', 'pBurst', 'xVar', 'genotype', 'grpVar', 'unitType');


tblVivo = mcu_tblVivo('basepaths', basepaths, 'presets', {'spkStates'}, ...
    'flgClean', true);