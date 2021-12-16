% example calls for frequently used functions 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general boolean arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
forceA = false;     % force analyze 
forceL = false;     % force save
saveFig = false;
graphics = false;
saveVar = false;

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
% preprocessing of raw files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = 'J:\Data\lh99\2021-12-15_21-02-03';
mapch = [26,27,28,30,2,3,31,29,4,5,6,7,8,9,10,11,12,13,14,15,1,16,17,...
    32,18,19,20,21,22,23,24,25,33,34,35,36,37];
rmvch = [];
% mapch = [1 : 21];
% rmvch = [17];

% tank to dat
store = 'Raw1';
blocks = [2 : 8];
chunksize = 300;
clip = cell(1, 1);
datInfo = tdt2dat('basepath', basepath, 'store', store, 'blocks',  blocks,...
    'chunksize', chunksize, 'mapch', mapch, 'rmvch', rmvch, 'clip', clip);

% open ephys to dat
exp = [1];
rec = cell(max(exp), 1);
datInfo = preprocOE('basepath', basepath, 'exp', exp, 'rec', rec,...
    'rmvch', rmvch, 'mapch', mapch, 'concat', true,...
    'nchans', length(mapch), 'fsIn', 20000);

% digital input from OE
clear recPath
recPath{1} = 'J:\Data\lh99\2021-12-15_09-01-22\RecordNode107\experiment1\recording1';
recPath{2} = 'J:\Data\lh99\2021-12-15_09-01-22\RecordNode107\experiment3\recording1';
recPath{3} = 'J:\Data\lh99\2021-12-15_09-01-22\RecordNode107\experiment3\recording2';
exPathNew = pwd;
getDinOE('basepath', recPath, 'newpath', exPathNew,...
    'concat', true, 'saveVar', true);

% pre-process dat (remove channels, reorder, etc.)
clip = [];
datInfo = preprocDat('basepath', basepath, 'fname', '', 'mapch', mapch,...
    'rmvch', rmvch, 'nchans', length(mapch), 'saveVar', true,...
    'chunksize', 1e7, 'precision', 'int16', 'bkup', true,...
    'clip', clip);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LFP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create lfp file
LFPfromDat('basepath', basepath, 'cf', 450, 'chunksize', 5e6,...
    'nchans', 31, 'fsOut', 1250,...
    'fsIn', 20000)    

% load lfp
lfp = getLFP('basepath', basepath, 'ch', [spkpgrp{:}], 'chavg', {},...
    'fs', 1250, 'interval', [0 inf], 'extension', 'lfp',...
    'savevar', true, 'forceL', true, 'basename', '');

% inter ictal spikes
binsize = (2 ^ nextpow2(30 * lfp.fs));
iis = getIIS('sig', double(lfp.data(:, 1)), 'fs', lfp.fs, 'basepath', basepath,...
    'graphics', true, 'saveVar', true, 'binsize', binsize,...
    'marg', 0.05, 'basename', '', 'thr', [5 0], 'smf', 7,...
    'saveFig', false, 'forceA', true, 'spkw', false, 'vis', true);

% burst suppression
vars = {'std', 'max', 'sum'};
bs = getBS('sig', double(lfp.data(:, ch)), 'fs', lfp.fs,...
    'basepath', basepath, 'graphics', true,...
    'saveVar', true, 'binsize', 1, 'BSRbinsize', binsize, 'smf', smf,...
    'clustmet', 'gmm', 'vars', vars, 'basename', '',...
    'saveFig', false, 'forceA', true, 'vis', true);

% create 'emg' signal from lfp cross correlation
emglfp = getEMGfromLFP(double(lfp.data(:, :)),...
    'emgFs', 10, 'saveVar', true);

% create emg signal from accelerometer data
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
ss = as_wrapper(EEG, EMG, sigInfo, 'basepath', pwd, 'calfile', [],...
    'viaGui', false, 'forceCalibrate', true, 'inspectLabels', false,...
    'saveVar', true, 'forceAnalyze', true, 'fs', 1250, 'netfile', netfile,...
    'graphics', true);

% inspect separation after classifying / manual scoring
as_stateSeparation(EEG, EMG, labels)

% calc psd in states
EEG = double(bz_LoadBinary([basename, '.emg.dat'], 'duration', Inf,...
    'frequency', 3051.7578125, 'nchannels', 2, 'start', 0,...
    'channels', 2, 'downsample', 1));
load([basename, '.AccuSleep_labels.mat'])
[psdStates, faxis, emgRMS] = psd_states('eeg', EEG, 'emg', [],...
    'labels', labels, 'fs', 3051.7578125, 'graphics', true, 'sstates', [1, 4]);

% get confusion matrix between two labels
[ss.netPrecision, ss.netRecall] = as_cm(labels1, labels2);
       
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
    'graphics', true, 'saveFig', false,...
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

% firing rate
load([basename, '.spikes.cellinfo.mat'])
binsize = 60;
winBL = [0 30 * 60];
fr = firingRate(spikes.times, 'basepath', basepath, 'graphics', true, 'saveFig', false,...
    'binsize', binsize, 'saveVar', true, 'smet', 'MA', 'winBL',...
    winBL, 'winCalc', [0, Inf]);

% plot fr vs. time
plot_FRtime_session('basepath', pwd, 'grp', [1 : 4, 7, 8],...
    'frBoundries', [0.1 Inf; 0.1 Inf], 'muFlag', false)

% cluster validation
spikes = cluVal('spikes', spikes, 'basepath', basepath, 'saveVar', true,...
    'saveFig', false, 'force', true, 'mu', [], 'graphics', false,...
    'vis', 'on', 'spkgrp', spkgrp);

% pyr vs. int
cc = cellclass('basepath', basepath,...
    'waves', cat(1, spikes.rawWaveform{:})', 'saveVar', true,...
    'graphics', false, 'fs', fs);

% spike timing metrics
st = spktimesMetrics('winCalc', ss.stateEpochs([1, 4]));

% organize firing rate 
[mfrCell, gainCell] = org_mfrCell('spikes', spikes, 'cm', cm, 'fr', fr,...
    'timebins', [1 Inf], 'dataType', 'su', 'grp', [1 : 4], 'suFlag', suFlag,...
    'frBoundries', [0 Inf; 0 Inf], 'stateIdx', stateIdx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% utilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert binary vector to epochs
epochs = binary2epochs('vec', binaryVec, 'minDur', 20, 'maxDur', 100,...
    'interDur', 30, 'exclude', false);

% convert cell to mat padded by nan (for different length arrays)
mat = cell2nanmat(c, dim);

% get size of bytes for a variable or class
nbytes = class2bytes(x, 'var', var);

% convert number to chunks. can clip parts or apply an overlap
chunks = n2chunks('n', nsamps, 'chunksize', chunksize, 'clip', clip);

% basically a wrapper for histcounts
[rate, binedges, tstamps, binidx] = times2rate(spktimes,...
    'binsize', binsize, 'winCalc', winCalc, 'c2r', true);

% convert basename to date time 
[~, dtFormat] = guessDateTime(dtstr);

% convert timestamp to datetime based on the date time string
[dt, tstamp] = tstamp2time(varargin);

% select specific units based on firing rate, cell class, etc.
units = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, unitClass);

% creates a cell array of variables
varArray = getSessionVars('dirnames', {basename}, 'mousepath', mousepath,...
    'sortDir', false);

% assigns vars to base workspace
assignVars(varArray, isession)

% set matlab graphics to custom or factory defaults
setMatlabGraphics(false)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis across sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fEPSP_sessions
fr_sessions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fEPSP signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see slutsky_fepsp repository 

% fEPSP from Dat
intens = [40 80 100 150 200];
fepsp = fEPSPfromDat('basepath', filepath, 'fname', '', 'nchans', nchans,...
    'spkgrp', spkgrp, 'intens', intens, 'concat', false, 'saveVar', false,...
    'force', true, 'extension', 'dat', 'recSystem', 'oe',...
    'protocol', 'io', 'graphics', true);

% fEPSP from WCP
basepath = 'C:\Users\heiml\Downloads\fEPSP\fEPSP\lh60';
intens = [30 : 10 : 60];
sfiles = [];
fepsp = fEPSPfromWCP('basepath', basepath, 'sfiles', [],...
    'sufx', 'io1', 'force', true, 'protocol', 'io',...
    'intens', intens, 'inspect', true, 'fs', 20000);

