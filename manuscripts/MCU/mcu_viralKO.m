



%% ========================================================================
%  CAG:MCU-KO - METRICS (st & bursts)
%  ========================================================================
% Three MCU-KO recordings (CAG cohort, 'ra' path key; TDT, fs ~24.4 kHz).
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
    
    % Sleep Signal
    sSig = as_prepSig([basename, '.lfp'], [basename, '.emg.dat'],...
        'eegCh', [1 : 4], 'emgCh', 1, 'saveVar', true, 'emgNchans', 4, 'eegNchans', 12,...
        'inspectSig', false, 'forceLoad', true, 'eegFs', 1250, 'emgFs', 3051.76,...
        'emgCf', [80 450]);
    
    labelsmanfile = [basename, '.sleep_labelsMan.mat'];
    AccuSleep_viewer(sSig, labels, labelsmanfile)
    
    % classify with a network
    netfile = 'D:\Code\slutsky_ECInVivo\lfp\SleepStates\AccuSleep\trainedNetworks\net_230212_103132.mat';
    calData = [];
    ss = as_classify(sSig, 'basepath', basepath, 'inspectLabels', false,...
        'saveVar', true, 'forceA', true, 'netfile', netfile,...
        'graphics', true, 'calData', calData);

    % spktimes
    spktimes = v(iFile).spikes.times;

    % Spike timing metrics
    st = spktimes_metrics('spktimes', spktimes, 'sunits', [], ...
        'bins', {[0, Inf]}, 'flgForce', true, 'flgSave', true, ...
        'flgAll', false);

    % Burst detection
    burst = burst_detect(spktimes, ...
        'minSpks', minSpks, ...
        'isiStart', isiStart, ...
        'isiEnd', isiEnd * 2, ...
        'minDur', minDur, ...
        'minIBI', minIBI, ...
        'flgForce', true, 'flgSave', true, 'flgPlot', false);

    % Burst statistics
    stats = burst_stats(burst, spktimes, 'winCalc', [], 'flgSave', true);

    % Ripples
    ripp = ripp_wrapper('basepath', pwd, ...
        'win', [0 12] * 3600, ...
        'rippCh', [1 : 4], ...
        'flgPlot', true, ...
        'flgSave', true, ...
        'flgNS', true, ...
        'flgForce', false);
    
    % Epileptiform discharges
    ed = ed_wrapper('basepath', basepath, ...
        'flgSave', true, ...
        'flgPlot', true, ...
        'flgForce', true);  


end


[cfgData, cfgGui] = guiPath_presets('ed');
cfgData = guiPath_load(cfgData);
[hFig, cfgData] = guiPath(basepath, 'cfgData', cfgData, 'cfgGui', cfgGui);




%% ========================================================================
%  CAG:MCU-KO - CLASSIFICATION
%  ========================================================================
% Assign RS / FS unit types and write units.mat. Two methods, toggled by flgCE:
%   flgCE = false : my pipeline. utypes_classify (GMM on waveform features)
%                   opens utypes_gui. Inspect the scatter / waveforms, reassign
%                   points, then click "Push Units" to save units.mat.
%   flgCE = true  : fall back to CellExplorer putativeCellType (Pyramidal -> RS,
%                   Narrow Interneuron -> FS, anything else -> Other).

flgCE = true;

basepaths = mcu_basepaths('ra');
nFiles = length(basepaths);

if ~flgCE

    % My pipeline: GMM classification, then manual curation in the GUI
    fetSelect = {'Asym', 'Hpk', 'TP'};
    tblUnit = utypes_classify('basepaths', basepaths, 'fetSelect', fetSelect, ...
        'rsPrior', 0.97, 'regVal', 0.01, 'flgPlot', true);

else

    % Fallback: CellExplorer putativeCellType
    v = basepaths2vars('basepaths', basepaths, 'vars', {'cell_metrics'});

    for iFile = 1 : nFiles

        % file
        basepath = basepaths{iFile};
        [~, basename] = fileparts(basepath);

        % Unit types from CellExplorer
        ctype = v(iFile).cell_metrics.putativeCellType(:);
        nUnits = numel(ctype);
        typeNum = zeros(nUnits, 1);                                  % Other
        typeNum(contains(ctype, 'Pyramidal'))          = 1;         % RS
        typeNum(contains(ctype, 'Narrow Interneuron')) = 2;         % FS

        units = struct();
        units.clean = false(2, nUnits);
        units.clean(1, :) = (typeNum == 1)';
        units.clean(2, :) = (typeNum == 2)';
        units.type = categorical(typeNum, [0, 1, 2], {'Other', 'RS', 'FS'});
        units.date = datetime('now');
        save(fullfile(basepath, [basename, '.units.mat']), 'units')

    end
end


%% ========================================================================
%  CAG:MCU-KO - BURSTINESS COMPARISON
%  ========================================================================
% Compare baseline burstiness of the three CAG:MCU-KO mice against the wt and
% mcu baseline groups. RS units only; CAG:MCU-KO is kept as its own group.

% CAG:MCU-KO group
basepaths = mcu_basepaths('ra');
tblCag = mcu_tblVivo('basepaths', basepaths, 'presets', {'burst'}, ...
    'flgClean', false);
tblCag.genotype = addcats(tblCag.genotype, {'CAG:MCU-KO'});
tblCag.genotype(:) = 'CAG:MCU-KO';

% Reference baseline groups (wt and mcu)
basepathsRef = [mcu_basepaths('wt_bsl'), mcu_basepaths('mcu_bsl')];
tblRef = mcu_tblVivo('basepaths', basepathsRef, 'presets', {'burst'}, ...
    'flgClean', false);
tblRef.genotype = addcats(tblRef.genotype, {'CAG:MCU-KO'});

% Combine
tblBrst = [tblRef; tblCag];
tblBrst.genotype = reordercats(tblBrst.genotype, ...
    {'Control', 'MCU-KO', 'CAG:MCU-KO'});

% Compare (switch yVar in the GUI for frBurst, br, bSize, etc.)
guiTbl_bar(tblBrst, 'yVar', 'pBurst', 'xVar', 'genotype', 'grpVar', 'unitType');







%% ========================================================================
%  TEST GUI
%  ========================================================================

basepaths = unique([mcu_basepaths('wt_bsl'), mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl')]);

basepath = basepaths{3};

% Signals (loaded once; sSig/specAdapter are full-session for the GUI)
[sig, emg, emgRms, fs, specAdapter, sSig] = ed_sigLoad(basepath);

ed = ed_wrapper('basepath', basepath, 'flgSave', true, 'flgPlot', false);  % save ed.mat for the EDs preset; suppress its auto-GUI

tic
guiPath(basepath, 'preset', 'EDs');
toc


guiPath(basepath);


tic
AccuSleep_viewer(sSig, [], [])
toc



