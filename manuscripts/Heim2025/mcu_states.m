
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pipeline for classification of a single mouse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create lfp file
LFPfromDat('basepath', basepath, 'cf', 450, 'chunksize', 5e6,...
    'nchans', nchans, 'fsOut', 1250, 'fsIn', fs)    


% create emg signal from accelerometer data
acc = EMGfromACC('basepath', basepath, 'fname', [basename, '.lfp'],...
    'nchans', nchans, 'ch', nchans - 2 : nchans, 'saveVar', true, 'fsIn', 1250,...
    'graphics', false, 'force', true);


% call for acceleration
sSig = as_prepSig([basename, '.lfp'], acc.mag,...
    'eegCh', [9 : 12], 'emgCh', [], 'saveVar', true, 'emgNchans', [],...
    'eegNchans', nchans, 'inspectSig', false, 'forceLoad', true,...
    'eegFs', 1250, 'emgFs', 1250, 'eegCf', [], 'emgCf', [10 450], 'fs', 1250);


% get spec outliers
otl = get_otlSpec('basepath', pwd, 'saveVar', true,...
    'flgForce', true, 'graphics', false, 'flgInspect', true);


% manually create labels. At this point, I manually classify at least 1
% hour of data for calibrating the signal and for future examination of the
% classification accuracy. 
labelsmanfile = [basename, '.sleep_labelsMan.mat'];
AccuSleep_viewer(sSig, labels, labelsmanfile)


% classify with a network
netfile = 'D:\Code\slutsky_ECInVivo\lfp\SleepStates\AccuSleep\trainedNetworks\net_230212_103132.mat';
ss = as_classify(sSig, 'basepath', basepath, 'inspectLabels', false,...
    'saveVar', true, 'forceA', true, 'netfile', netfile,...
    'graphics', true, 'calData', []);


% create bouts
minDur = 5;
interDur = 3;
bouts = as_bouts('labels', ss.labels,...
    'minDur', minDur, 'interDur', interDur, 'flgOtl', true,...
    'sstates', 1 : cfg.nstates, 'graphics', false);


% calculate EMG state bouts
minDur = 5;
interDur = 3;
ssEmg = as_emg('basepath', basepath, 'flgInspct', false,...
    'interDur', interDur, 'minDur', minDur);















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manually curate state bouts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

queryStr = 'wt_bsl';
basepaths = [mcu_sessions(queryStr)];
nfiles = length(basepaths);

% select session
ifile = 3;
basepath = basepaths{ifile};
[~, basename] = fileparts(basepath);
cd(basepath)

% load manual labels
if exist(fname_labelsMan, "file")

    % create bkup of manual labels
    mkdir(fullfile(basepath, 'ss'))
    fname_labelsMan = [basename, '.sleep_labelsMan.mat'];
    ctime = datetime('now', 'Format', 'yyMMdd_HHmmss');
    bkupName = [basename, '.labelsMan', '_', char(ctime), '.mat'];
    fname_bkup = fullfile(basepath, 'ss', bkupName);
    copyfile(fname_labelsMan, fname_bkup)

else
    load([basename, '.sleep_states.mat'])
    labels = ss.labels;
end

% load sleep sig
sSig = load([basename, '.sleep_sig.mat']);

% manually revise network labels
AccuSleep_viewer(sSig, labels, fname_labelsMan)








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reclassify states automatically
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepaths = [mcu_sessions('wt_wsh'), mcu_sessions('mcu_wsh')];
% basepaths = [mcu_sessions('lh100')];
nfiles = length(basepaths);

for ifile = 1 : nfiles
    basepath = basepaths{ifile};
    [~, basename] = fileparts(basepath);
    cd(basepath)

    % backup labels man file
    fname_labelsMan = [basename, '.sleep_labelsMan.mat'];
    ctime = datetime('now', 'Format', 'yyMMdd_HHmmss');
    bkupName = [basename, '.labelsMan', '_', char(ctime), '.mat'];
    fname_bkup = fullfile(basepath, 'ss', bkupName);
    if exist("fname_labelsMan", 'file')
        copyfile(fname_labelsMan, fname_bkup)
    end

    % backup labels rev file and rename to man
    fname_labelsRev = fullfile(basepath, [basename, '.sleep_labelsRev.mat']);
    if exist(fname_labelsRev, 'file')
        bkupName = [basename, '.labelsRev', '_', char(ctime), '.mat'];
        fname_bkup = fullfile(basepath, 'ss', bkupName);
        copyfile(fname_labelsRev, fname_bkup)
        movefile(fname_labelsRev, fname_labelsMan);
    end

    % load data
    ssigfile = fullfile(basepath, [basename, '.sleep_sig.mat']);
    sSig = load(ssigfile);

    % classify with a network
    netfile = 'D:\Code\slutsky_ECInVivo\lfp\SleepStates\AccuSleep\trainedNetworks\net_230212_103132.mat';
    ss = as_classify(sSig, 'basepath', basepath, 'inspectLabels', false,...
        'saveVar', true, 'forceA', true, 'netfile', netfile,...
        'graphics', true, 'calData', []);

    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find spec outliers in each file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mice = [mcu_sessions('mcu')];
for imouse = 1 : length(mice)

    mname = mice{imouse};
    basepaths = mcu_sessions(mname);
    nfiles = length(basepaths);

    for ifile = 1 : nfiles
        basepath = basepaths{ifile};
        cd(basepath)
        otl = get_otlSpec('basepath', pwd, 'saveVar', true,...
            'flgForce', true, 'graphics', false, 'flgInspect', true);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recalculate state bouts (accusleep)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mice = [mcu_sessions('mcu')];
vars = ["sleep_states"];
cfg = as_loadConfig();
minDur = 5;
interDur = 3;

for imouse = 1 : length(mice)

    mname = mice{imouse};
    [basepaths, v] = [mcu_sessions(mname, vars)];
    nfiles = length(basepaths);
    mpath = fileparts(basepaths{1});

    for ifile = 1 : nfiles

        % files
        basepath = basepaths{ifile};
        [~, basename] = fileparts(basepath);
        cd(basepath)
        fname_sig = fullfile(basepath, [basename, '.sleep_sig.mat']);
        fname_states = [basename '.sleep_states.mat'];

        % create bouts
        ss = v(ifile).ss;
        bouts = as_bouts('labels', ss.labels,...
            'minDur', minDur, 'interDur', interDur, 'flgOtl', true,...
            'sstates', 1 : cfg.nstates, 'graphics', false);

        % inspect new labels
        flg_inspect = false;
        if flg_inspect
            sSig = load(fname_sig);
            AccuSleep_viewer(sSig, bouts.labels, [])
            AccuSleep_viewer(sSig, bouts.labelsOrig, [])
        end

        % update ss.bouts and save
        flg_update = true;
        if flg_update
            ss.bouts = bouts;
            save(fname_states, 'ss');
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate state bouts (emg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mice = [mcu_sessions('mcu')];
cfg = as_loadConfig('flgEmg', true);
minDur = 5;
interDur = 3;

for imouse = 1 : length(mice)

    mname = mice{imouse};
    [basepaths] = [mcu_sessions(mname)];
    nfiles = length(basepaths);
    mpath = fileparts(basepaths{1});

    for ifile = 1 : nfiles

        % files
        basepath = basepaths{ifile};
        [~, basename] = fileparts(basepath);
        cd(basepath)
        fname_sig = fullfile(basepath, [basename, '.sleep_sig.mat']);
        fname_states = [basename '.sleep_statesEmg.mat'];

        % create bouts
        ssEmg = as_emg('basepath', basepath, 'flgInspct', false,...
            'interDur', interDur, 'minDur', minDur);

        % inspect new labels
        flg_inspect = false;
        if flg_inspect
            sSig = load(fname_sig);
            AccuSleep_viewer(sSig, ssEmg.bouts.labels, [], cfg)
            AccuSleep_viewer(sSig, ssEmg.bouts.labelsOrig, [], cfg)
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LME 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data for each group
grps = {'wt', 'mcu'};
idx_rmDays = [];              % remove bac on and off
clear grppaths
for igrp = 1 : length(grps)

    mnames = mcu_sessions(grps(igrp));  
    for imouse = 1 : length(mnames)
        basepaths = mcu_sessions(mnames{imouse});
        basepaths(idx_rmDays) = [];
        grppaths{igrp}(imouse, :) = string(basepaths)';
    end
end


% Bout length of WT vs. MCU across time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frml = 'BDur ~ Group * Day + (1|Mouse)';

% organize for lme
var_field = '';
[lme_tbl, lme_cfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flg_emg', true, 'var_field', var_field, 'vCell', {});

% normalize to baseline
plot_tbl = lme_normalize('lme_tbl', lme_tbl, 'normVar', 'Day', 'groupVars', {'Group', 'State'});

% select state
istate = categorical({'Low EMG'});
plot_tbl = plot_tbl(plot_tbl.State == istate, :);

% run lme
lme = fitlme(plot_tbl, lme_cfg.frml);

% plot
fh = lme_plot(plot_tbl, lme, 'ptype', 'line');





% save
frml = [char(lme.Formula), '_', char(istate)];
frml = frml2char(frml, 'rm_rnd', false);
th = get(fh, 'Children'); % if you're sure it's a tiled layout
title(th, frml, 'interpreter', 'none')
grph_save('fh', fh, 'fname', frml, 'frmt', {'ai', 'jpg'})

[prism_data] = fh2prism(fh);

% test interaction of wt and mcu across all days
interactionIndices = contains(lme.CoefficientNames, 'Group_MCU-KO:Day'); % Logical vector for interaction terms
contrastVector = zeros(1, length(lme.Coefficients.Estimate));
contrastVector(interactionIndices) = 1; % Test only interaction terms

% Define contrast vector for the difference between MCU-KO and WT during WASH
contrastVector = zeros(1, length(lme.Coefficients.Estimate));
contrastVector(strcmp(lme.CoefficientNames, 'Group_MCU-KO')) = 1; % Main effect of Group
contrastVector(strcmp(lme.CoefficientNames, 'Group_MCU-KO:Day_WASH')) = 1; % Interaction term

% Run the test
[p, h, stat] = coefTest(lme, contrastVector)

% Bout length of High vs. Low EMG across time per group
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% organize for lme
frml = 'BLen ~ Day + (1|Mouse)';
[lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, true);

% select group
igrp = categorical({'WT'});
plot_tbl = lme_tbl(lme_tbl.Group == igrp, :);

% select state
istate = categorical({'Low EMG'});
plot_tbl = lme_tbl(plot_tbl.State == istate, :);

% run lme
lme = fitlme(plot_tbl, lme_cfg.frml);

% plot
fh = mcu_lmePlot(plot_tbl, lme, 'ptype', 'line');

% Bout length of WT vs. MCU across states during baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grps = {'wt_bsl'; 'mcu_bsl'};
grps = {'wt_wsh'; 'mcu_wsh'};

clear grppaths
for igrp = 1 : length(grps)
    grppaths{igrp} = string(mcu_sessions(grps{igrp})');
end

% organize for lme
frml = 'BLen ~ Group * State + (1|Mouse)';
[lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, false);

% run lme
plot_tbl = lme_tbl;
lme = fitlme(plot_tbl, lme_cfg.frml);

% plot
fh = mcu_lmePlot(plot_tbl, lme, 'ptype', 'line');

% save
grph_save('fh', fh, 'fname', [char(lme.Formula), '_AS'], 'frmt', {'ai', 'jpg'})

[prism_data] = fh2prism(fh);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare mcu and wt during baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% can use states either from accusleep or emg 
flg_emg = true;

% can use bouts after cleaning and merging or original
flg_clean = true;

% organize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data
clear basepaths
basepaths{1} = [mcu_sessions('wt_bsl')];
basepaths{2} = [mcu_sessions('mcu_bsl')];
basepaths{2}(end) = [];
nfiles = max(cellfun(@length, basepaths));

if flg_emg
    vars = ["sleep_statesEmg"];
    sstates = [1, 2];
else
    vars = ["sleep_states"];
    sstates = [1, 4, 5];
end
cfg = as_loadConfig('flgEmg', flg_emg);
nstates = length(sstates);
clr = cfg.colors(sstates);
snames = cfg.names(sstates);

% initialize
boutLen = nan(2, nfiles, length(sstates));
prctDur = nan(2, nfiles, length(sstates));

for igrp = 1 : 2
    v = basepaths2vars('basepaths', basepaths{igrp}, 'vars', vars);
    nfiles = length(basepaths{igrp});
    ss_sfx = fieldnames(v);

    ss = catfields([v(:).(ss_sfx{1})], 1, true);   
    for ifile = 1 : nfiles        
        
        if flg_clean
            btimes = ss.bouts.times(ifile, sstates);
        else
            btimes = ss.bouts.timesOrig(ifile, sstates);
        end
        winlen = diff(ss.bouts.info.timewins(ifile, :));

        % bout length
        boutLen_cell = cellfun(@(x) (diff(x')'), btimes, 'uni', false);
        boutLen(igrp, ifile, :) = mean(cell2padmat(boutLen_cell, 2), 'omitnan');

        % state duration (%)
        totDur = cellfun(@sum, boutLen_cell);
        prctDur(igrp, ifile, :) = totDur * 100 ./ winlen;
    end
end


% Plot Comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open figure
fh = figure;
setMatlabGraphics(true);
set(fh, 'WindowState','maximized');
tlayout = [2, length(sstates)];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';

% Plot bout length for each state
tilebias = 0;
for istate = 1 : length(sstates)
    axh = nexttile(th, tilebias + istate); hold on; cla;
    stateIdx = sstates(istate);
    dataMat = boutLen(:, :, istate)';
    plot_boxMean('axh', axh, 'dataMat', dataMat, 'allPnt', false,...
        'plotType', 'bar');
    ylabel('Bout Length (s)');
    title([snames{istate}]);
    ylim([0, ceil(max(dataMat(:)))]);
    xlim([0.5 2.5])
    xticks([1, 2])
    xticklabels({'WT', 'MCU-KO'});
end


% Plot bout length for each state
tilebias = nstates;
for istate = 1 : length(sstates)
    axh = nexttile(th, tilebias + istate); hold on; cla;
    stateIdx = sstates(istate);
    dataMat = prctDur(:, :, istate)';
    plot_boxMean('axh', axh, 'dataMat', dataMat, 'allPnt', false,...
        'plotType', 'bar');
    ylabel('State Duration (%)');
    title([snames{istate}]);
    ylim([0, ceil(max(dataMat(:)))]);
    xlim([0.5 2.5])
    xticks([1, 2])
    xticklabels({'WT', 'MCU-KO'});
end
