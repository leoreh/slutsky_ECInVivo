 % preproc_wrapper
basepath = 'D:\LFP\SWC\TDT\TG\TG2\tg2_210716_2020';
cd(basepath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open ephys
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = 'J:\Data\lh99\2021-12-13_09-13-55';
rmvch = [];
mapch = [26,27,28,30,2,3,31,29,4,5,6,7,8,9,10,11,12,13,14,15,1,16,17,32,18,19,20,21,22,23,24,25,33,34,35,36,37];
exp = [1];
rec = cell(max(exp), 1);
% rec{1} = [2, 3];
datInfo = preprocOE('basepath', basepath, 'exp', exp, 'rec', rec,...
    'rmvch', rmvch, 'mapch', mapch, 'concat', true,...
    'nchans', length(mapch), 'fsIn', 20000);

% pre-process dat
datInfo = preprocDat('basepath', 'J:\Data\lh99\lh99_211211_090453',...
    'fname', '', 'mapch', mapch,...
    'rmvch', rmvch, 'nchans', nchans, 'saveVar', true,...
    'chunksize', 5e6, 'precision', 'int16', 'bkup', false);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tdt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = 'J:\Data\lh96_211208_070600';
store = 'Raw2';
blocks = [2, 6, 7];
chunksize = 300;
mapch = [1 : 4];
mapch = [1 : 16];
rmvch = [2, 4];
rmvch = [7];
clip = cell(1, 1);
datInfo = tdt2dat('basepath', basepath, 'store', store, 'blocks',  blocks,...
    'chunksize', chunksize, 'mapch', mapch, 'rmvch', rmvch, 'clip', clip);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% session info (cell explorer foramt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'force', true, 'saveVar', true);      
basepath = session.general.basePath;
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;
[~, basename] = fileparts(basepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spike sorting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in this pipeline waveforms are extracted from the whitened data and are
% used all the way until manual curation (including calculation of
% isolation distance). Afterwards, cell explorer resnips the waveform for
% the dat file and detrends them such that cell class is determined
% from the raw waveforms.

% spike detection from temp_wh
[spktimes, ~] = spktimesWh('basepath', basepath, 'fs', fs, 'nchans', nchans,...
    'spkgrp', spkgrp, 'saveVar', true, 'saveWh', true,...
    'graphics', false, 'force', true);

% spike rate per tetrode. note that using firingRate requires
% special care becasue spktimes is given in samples and not seconds
for igrp = 1 : length(spkgrp)
    spktimes{igrp} = spktimes{igrp} / fs;
end
sr = firingRate(spktimes, 'basepath', basepath,...
    'graphics', false, 'saveFig', false,...
    'binsize', 60, 'saveVar', 'sr', 'smet', 'none',...
    'winBL', [0 Inf]);

% create ns files 
dur = [];
t = [];
spktimes2ns('basepath', basepath, 'fs', fs,...
    'nchans', nchans, 'spkgrp', spkgrp, 'mkClu', true,...
    'dur', dur, 't', t, 'psamp', [], 'grps', [1 : length(spkgrp)],...
    'spkFile', 'temp_wh');

% clip ns files
nsClip('dur', -420, 't', [], 'bkup', true, 'grp', [3 : 4]);

% clean clusters after sorting 
cleanCluByFet('basepath', pwd, 'manCur', true, 'grp', [1 : 4])

% cut spk from dat and realign
fixSpkAndRes('grp', [1 : 4], 'dt', 0, 'stdFactor', 0, 'resnip', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spikes (post-sorting)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cell explorer metrics
cell_metrics = ProcessCellMetrics('session', session,...
    'manualAdjustMonoSyn', false, 'summaryFigures', false,...
    'debugMode', true, 'transferFilesFromClusterpath', false,...
    'submitToDatabase', false, 'getWaveformsFromDat', true,...
    'forceReloadSpikes', false);
% cell_metrics = CellExplorer('basepath', pwd);

% load spikes ce format
% spikes = loadSpikes('format', 'klustakwik', 'getWaveformsFromSource', true, 'LSB', 1);

% cluster validation
load([basename, '.spikes.cellinfo.mat'])
spikes = cluVal('spikes', spikes, 'basepath', basepath, 'saveVar', true,...
    'saveFig', false, 'force', true, 'mu', [], 'graphics', false,...
    'vis', 'on', 'spkgrp', spkgrp);

% pyr vs. int
cc = cellclass('basepath', basepath,...
    'waves', cat(1, spikes.rawWaveform{:})', 'saveVar', true,...
    'graphics', false, 'fs', fs);

% spike timing metrics
st = spktimesMetrics('winCalc', ss.boutTimes([1, 4]));

% firing rate
binsize = 60;
winBL = [0 30 * 60];
fr = firingRate(spikes.times, 'basepath', basepath, 'graphics', false, 'saveFig', false,...
    'binsize', binsize, 'saveVar', true, 'smet', 'MA', 'winBL',...
    winBL, 'winCalc', [0, Inf]);

[mfrCell, gainCell] = org_mfrCell('spikes', spikes, 'cm', cm, 'fr', fr,...
    'timebins', [1 Inf], 'dataType', 'su', 'grp', [1 : 4], 'suFlag', suFlag,...
    'frBoundries', [0 Inf; 0 Inf], 'stateIdx', stateIdx);

% CCG
binSize = 0.001; dur = 0.12; % low res
binSize = 0.0001; dur = 0.02; % high res
[ccg, t] = CCG({xx.times{:}}, [], 'duration', dur, 'binSize', binSize, 'norm', 'counts');
u = 20;
plotCCG('ccg', ccg(:, u, u), 't', t, 'basepath', basepath,...
    'saveFig', false, 'c', {'k'}, 'u', spikes.UID(u));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lfp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create lfp
LFPfromDat('basepath', basepath, 'cf', 450, 'chunksize', 5e6,...
    'nchans', nchans, 'fsOut', 1250,...
    'fsIn', fs)   

% load lfp
lfp = getLFP('basepath', basepath, 'ch', [spkgrp{:}], 'chavg', {},...
    'fs', 1250, 'interval', [0, Inf], 'extension', 'lfp',...
    'savevar', true, 'forceL', true, 'basename', '');

% remove 50 Hz from signal
emgOrig = filterLFP(emgOrig, 'fs', 1250, 'stopband', [49 51],...
    'dataOnly', true, 'saveVar', false, 'graphics', false);

% create emg signal from accelerometer
acc = EMGfromACC('basepath', basepath, 'fname', [basename, '.lfp'],...
    'nchans', nchans, 'ch', nchans - 2 : nchans, 'saveVar', true, 'fsIn', 1250,...
    'graphics', false);

% get ripples
ripp = getRipples('basepath', basepath, 'rippCh', [23], 'emgCh', [33],...
    'emg', [], 'recWin', [0, 4 * 60 * 60]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sleep states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prep signal
[EMG, EEG, sigInfo] = as_prepSig([basename, '.lfp'], [basename, '.lfp'],...
    'eegCh', [25 : 28], 'emgCh', [33], 'saveVar', true, 'emgNchans', [],...
    'eegNchans', 37, 'inspectSig', false, 'forceLoad', true,...
    'eegFs', 1250, 'emgFs', 1250, 'emgCf', [10 200]);

% manually create labels
labelsmanfile = [basename, '.AccuSleep_labelsMan.mat'];
AccuSleep_viewer(EEG, EMG,  1250, 1, [], labelsmanfile)

% classify with a network
netfile = 'D:\Code\slutskycode\extracellular in vivo\lfp\SleepStates\AccuSleep\trainedNetworks\net_211107_141044.mat';
% netfile = [];
ss = as_wrapper(EEG, EMG, sigInfo, 'basepath', pwd, 'calfile', [],...
    'viaGui', false, 'forceCalibrate', true, 'inspectLabels', false,...
    'saveVar', true, 'forceAnalyze', true, 'fs', 1250, 'netfile', netfile,...
    'graphics', true);

% inspect separation after classifying / manual scoring
as_stateSeparation(EEG, EMG, labels)

% calc psd in states
EEG = double(bz_LoadBinary([basename, '.lfp'], 'duration', Inf,...
    'frequency', 1250, 'nchannels', nchans, 'start', 0,...
    'channels', [4 : 7], 'downsample', 1));
EEG = mean(EEG, 2);

EEG = double(bz_LoadBinary([basename, '.emg.dat'], 'duration', Inf,...
    'frequency', 3051.7578125, 'nchannels', 2, 'start', 0,...
    'channels', 2, 'downsample', 1));

load([basename, '.AccuSleep_labels.mat'])
[psdStates, faxis, emgRMS] = psd_states('eeg', EEG, 'emg', [],...
    'labels', labels, 'fs', 3051.7578125, 'graphics', true, 'sstates', [1, 4]);

% get confusion matrix between two labels
[ss.netPrecision, ss.netRecall] = as_cm(labels1, labels2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% handle dat  files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% cat dat
nchans = 37;
fs = 20000;
newpath = 'J:\Data\lh98';
datFiles{1} = 'J:\Data\lh98\lh98_211219_085802\lh98_211219_085802.dat';
datFiles{2} = 'J:\Data\lh98\lh98_211219_181302\lh98_211219_181302.dat';
sigInfo = dir(datFiles{1});
nsamps = floor(sigInfo.bytes / class2bytes('int16') / nchans);
parts{1} = round([0, Inf]);
parts{2} = round([0,  (748 * 60 + 8) * fs]);

catDatMemmap('datFiles', datFiles, 'newpath', newpath, 'parts', [],...
    'nchans', nchans, 'saveVar', true)

clear source
source{1} = 'J:\Data\lh98\lh98_211220_104619\lh98_211220_104619.dat';
source{2} = 'J:\Data\lh98\lh98_211220_211037\lh98_211220_211037.dat';
destination = 'J:\Data\lh98\lh98_211220.dat';
cmd = ['!copy /b ' strjoin(source, ' + '), ' ' destination];
eval(cmd);

% preproc dat
clip = [(1290 * 60 + 17) * fs, Inf];
datInfo = preprocDat('basepath', pwd,...
    'fname', [basename, '.dat'], 'mapch', 1 : nchans,...
    'rmvch', [], 'nchans', nchans, 'saveVar', false, 'clip', clip,...
    'chunksize', 5e7, 'precision', 'int16', 'bkup', true);


