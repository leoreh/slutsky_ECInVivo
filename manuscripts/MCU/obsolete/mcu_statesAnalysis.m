


%% ========================================================================
%  ANALYSIS PIPELINE (per Mouse)
%  ========================================================================

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



%% ========================================================================
%  MANUAL CURATION
%  ========================================================================

queryStr = 'wt_bsl';
basepaths = mcu_basepaths(queryStr);
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





%% ========================================================================
%  AUTOMATIC RE-CLASSIFICATION
%  ========================================================================

basepaths = [mcu_basepaths('wt_wsh'); mcu_basepaths('mcu_wsh')];
% basepaths = mcu_basepaths('lh100');
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


%% ========================================================================
%  SPEC OUTLIERS
%  ========================================================================

mice = unique(get_mname(mcu_basepaths('mcu')));
for imouse = 1 : length(mice)

    mname = mice{imouse};
    basepaths = mcu_basepaths(mname);
    nfiles = length(basepaths);

    for ifile = 1 : nfiles
        basepath = basepaths{ifile};
        cd(basepath)
        otl = get_otlSpec('basepath', pwd, 'saveVar', true,...
            'flgForce', true, 'graphics', false, 'flgInspect', true);
    end
end


%% ========================================================================
%  STATE BOUTS (ACCUSLEEP)
%  ========================================================================

mice = unique(get_mname(mcu_basepaths('mcu')));
vars = ["sleep_states"];
cfg = as_loadConfig();
minDur = 5;
interDur = 3;

for imouse = 1 : length(mice)

    mname = mice{imouse};
    basepaths = mcu_basepaths(mname);
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);
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


%% ========================================================================
%  STATE BOUTS (EMG)
%  ========================================================================

mice = unique(get_mname(mcu_basepaths('mcu')));
cfg = as_loadConfig('flgEmg', true);
minDur = 5;
interDur = 3;

for imouse = 1 : length(mice)

    mname = mice{imouse};
    basepaths = mcu_basepaths(mname);
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



