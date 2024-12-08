
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manually curate state epochs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

queryStr = 'mcu_bsl';
basepaths = [mcu_sessions(queryStr)];
nfiles = length(basepaths);

% select session
ifile = 2;
basepath = basepaths{ifile};
[~, basename] = fileparts(basepath);
cd(basepath)

% create bkup of manual labels
labelsrevfile = [basename, '.sleep_labelsRev.mat'];
if exist(labelsrevfile, "file")
    load(labelsrevfile)
    mkdir(fullfile(basepath, 'ss'))
    ctime = datetime('now', 'Format', 'yyMMdd_HHmmss');
    bkupName = [basename, '.labelsRev', '_', char(ctime), '.mat'];
    bkupFile = fullfile(basepath, 'ss', bkupName);
    save(bkupFile, 'labels')
else
    load([basename, '.sleep_states.mat'])
    labels = ss.labels;
end

% load sleep sig 
sSig = load([basename, '.sleep_sig.mat']);

% manually revise network labels
AccuSleep_viewer(sSig, labels, labelsrevfile)

% create bkup of labelsMan
labelsmanfile = [basename, '.sleep_labelsMan.mat'];
ctime = datetime('now', 'Format', 'yyMMdd_HHmmss');
bkupName = [basename, '.labelsMan', '_', char(ctime), '.mat'];
bkupFile = fullfile(basepath, 'ss', bkupName);
copyfile(labelsmanfile, bkupFile)

% replace labelsMan with revised labels
idx = 1 : 6 * 60 * 60;
load(labelsrevfile)
labelsRev = labels;
labels = ones(length(labelsRev), 1) * 8;
labels(idx) = labelsRev(idx);
AccuSleep_viewer(sSig, labels, labelsmanfile)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare manual and automatic classification (ra data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

queryStr = 'mcu_bsl';
vars = ["sleep_states"];
[basepaths, v] = mcu_sessions(queryStr, vars);
nfiles = length(basepaths);


% params
interDur = 0;
minDur = 0;
sstates = [1, 4, 5];
nstates = length(sstates);

% file
ifile = 2;
basepath = basepaths{ifile};
cd(basepath)
[~, basename] = fileparts(basepath);
ss = v(ifile).ss;

idx = 1 : 6 * 60 * 60;
idx = 1 : length(ss.labels);
[stateEpochs, epochStats] = as_epochs('labels', ss.labels_man(idx),...
    'minDur', minDur, 'interDur', interDur);
[stateEpochs2, epochStats2] = as_epochs('labels', ss.labels_net(idx),...
    'minDur', minDur, 'interDur', interDur);

epLen = cell(nstates, 1);
for istate = 1 : nstates
    tmp{1} = epochStats.epLen{sstates(istate)};
    tmp{2} = epochStats2.epLen{sstates(istate)};
    epLen{istate} = cell2nanmat(tmp, 2);
end

prctDur = epochStats.prctDur(sstates);
prctDur(2, :) = epochStats2.prctDur(sstates);


% -------------------------------------------------------------------------
% compare state dependent firing rate
load([basename, '.spikes.cellinfo.mat'])
load([basename, '.units.mat'])

fr = calc_fr(spikes.times, 'basepath', basepath,...
    'graphics', false, 'binsize', 60, 'saveVar', false, 'forceA', true,...
    'smet', 'none', 'winBL', [0, Inf], 'winCalc', [0, Inf], 'stateEpochs', stateEpochs);

fr2 = calc_fr(spikes.times, 'basepath', basepath,...
    'graphics', false, 'binsize', 60, 'saveVar', false, 'forceA', true,...
    'smet', 'none', 'winBL', [0, Inf], 'winCalc', [0, Inf], 'stateEpochs', stateEpochs2);

iunit = 1;
unitIdx = units.clean(iunit, :);
clear mfr
for istate = 1 : nstates
    tmp{1} = fr.states.mfr(unitIdx, sstates(istate));
    tmp{2} = fr2.states.mfr(unitIdx, sstates(istate));
    mfr{istate} = cell2nanmat(tmp, 2);
end

% calculate % change
clear prctdiff
for istate = 1 : nstates
    prctdiff(:, istate) = (mfr{istate}(:, 1) - mfr{istate}(:, 2)) ./...
        ((mfr{istate}(:, 1) + mfr{istate}(:, 2)) / 2) * 100;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state (emg) epochs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for each session (and mouse), calculate EMG states
% calculate fr in timebins

% gen parans
saveVar = false;
graphics = true;

% state params
minDur = 10;
interDur = 3;
ftarget = [0.5 : 0.5 : 100];

% files
mname = 'lh134';
basepaths = [mcu_sessions(mname)];
nfiles = length(basepaths);

for ifile = 1 : nfiles

    basepath = basepaths{ifile};
    [~, basename] = fileparts(basepath);
    cd(basepath)


    [stateEpochs, epochStats] = as_epochs('minDur', minDur, 'interDur',...
        interDur, 'flgEmg', true, 'graphics', true, 'nbins', 2);
    
    % calc psd according to emg state separation
    psd = psd_states('basepath', basepath, 'sstates', [1, 2],...
        'sig', [], 'fs', fs, 'saveVar', saveVar,...
        'graphics', true, 'forceA', true, 'ftarget', ftarget,...
        'emgThr', [], 'flgEmg', true, 'stateEpochs', stateEpochs);

    fh = figure;
    plot_hypnogram('stateEpochs', psd.epochs.stateEpochs, 'clr', psd.info.clr, 'axh', gca,...
    'sstates', [1 : 2])

end
