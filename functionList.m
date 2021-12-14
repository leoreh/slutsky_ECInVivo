% list of functions and their calls

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general boolean arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
forceA = false;
forceL = false;
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
% IO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = 'D:\VMs\shared\lh57-lh69\lh57-lh69_200904';
store = 'Raw1';
fs = 24414.06;
blocks = [2 : 8];
chunksize = 300;
mapch = [1 : 16];
% mapch = [1 : 2 : 7, 2 : 2 : 8, 9 : 2 : 15, 10 : 2 : 16];
rmvch = [4];
clip = cell(1, 1);

% tank to dat
datInfo = tdt2dat('basepath', basepath, 'store', store, 'blocks',  blocks,...
    'chunksize', chunksize, 'mapch', mapch, 'rmvch', rmvch, 'clip', clip);

% open ephys to dat
exp = [10];
rec = cell(max(exp), 1);
datInfo = preprocOE('basepath', basepath, 'exp', exp, 'rec', rec,...
    'rmvch', rmvch, 'mapch', mapch, 'concat', true, 'nchans', 35,...
    'intens', intens);

% digital input from OE
getDinOE('basepath', recPath, 'newpath', exPathNew,...
    'concat', true, 'nchans', nchans, 'precision', 'int16',...
    'saveVar', true);

% pre-process dat (remove channels, reorder, etc.)
datInfo = preprocDat('basepath', basepath, 'fname', '', 'mapch', mapch,...
    'rmvch', rmvch, 'nchans', nchans, 'saveVar', true,...
    'chunksize', 1e7, 'precision', 'int16', 'bkup', true,...
    'clip', [808593667 Inf]);

% acceleration
newch = length(mapch) - length(rmvch);
chAcc = [newch : -1 : newch - 2];
EMGfromACC('basepath', exPathNew, 'fname', '',...
    'nchans', newch, 'ch', chAcc, 'force', false, 'saveVar', true,...
    'graphics', false, 'fsOut', 1250);

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

% sleep states ------------------------------------------------------------

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
st = spktimesMetrics('winCalc', ss.stateEpochs([1, 4]));

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
% utilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert binary vector to epochs
epochs = binary2epochs('vec', binaryVec, 'minDur', 20, 'maxDur', 100,...
    'interDur', 30, 'exclude', false);

% convert cell to mat padded by nan (for different length arrays)
mat = cell2nanmat(c, dim);

% get size of bytes for a variable or class
nbytes = class2bytes(x, 'var', var);

% convert basename to date time 
[dt, dtFormat] = guessDateTime(dtstr);

% convert number to chunks. can clip parts or apply an overlap
chunks = n2chunks('n', nsamps, 'chunksize', chunksize, 'clip', clip);

% basically a wrapper for histcounts
[rate, binedges, tstamps, binidx] = times2rate(spktimes,...
    'binsize', binsize, 'winCalc', winCalc, 'c2r', true);
        
% convert timestamp to datetime based on the date time string
[dt, tstamp] = tstamp2time(varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis across sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fEPSP_sessions
fr_sessions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fEPSP signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all functions here are obsolete. see slutsky_fepsp repository

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

% anesthesia states (see also aneStates_wrp)
[bs, iis, ep] = aneStates('ch', 1, 'basepath', basepath,...
    'basename', basename, 'graphics', graphics,...
    'saveVar', saveVar, 'saveFig', saveFig, 'forceA', forceA,...
    'binsize', 30, 'smf', 7, 'thrMet', 1);

