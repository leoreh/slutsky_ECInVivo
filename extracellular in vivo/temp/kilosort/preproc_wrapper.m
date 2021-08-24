 % preproc_wrapper
basepath = 'D:\VMs\shared\lh58\lh58_200901_080917';
cd(basepath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open ephys
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = 'D:\Data\lh93\2021-08-19_22-16-08';
rmvch = [4];
mapch = [1 : 19];
exp = [1];
rec = cell(max(exp), 1);
% rec{1} = [2, 3];
datInfo = preprocOE('basepath', basepath, 'exp', exp, 'rec', rec,...
    'rmvch', rmvch, 'mapch', mapch, 'concat', true,...
    'nchans', length(mapch), 'fsIn', 20000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tdt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = 'E:\Leore\LH\DataTanks\lh94\lh94_210823_121100';
store = 'Raw2';
blocks = [3];
chunksize = 300;
mapch = [1 : 4];
mapch = [1 : 16];
rmvch = [3, 4];
rmvch = [];
clip = cell(1, 1);
clip{1} = [0 seconds(hours(2) + minutes(41) + seconds(37))];
datInfo = tdt2dat('basepath', basepath, 'store', store, 'blocks',  blocks,...
    'chunksize', chunksize, 'mapch', mapch, 'rmvch', rmvch, 'clip', clip);

% get start/stop time of specific block from tdt
for iblock = 1 : length(blocks)
    blockpath = fullfile(basepath, ['block-', num2str(blocks(iblock))]);
    heads{iblock} = TDTbin2mat(blockpath, 'HEADERS', 1);
    heads{iblock}.startDate = datestr(datenum([1970, 1, 1, 0, 0, heads{iblock}.startTime]) + hours(3));
    heads{iblock}.stopDate = datestr(datenum([1970, 1, 1, 0, 0, heads{iblock}.stopTime]) + hours(3));
end

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
% field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intens = flip([100 : 20 : 200, 240, 280, 320]);
fepsp = fEPSPfromDat('basepath', basepath, 'fname', '', 'nchans', nchans,...
    'spkgrp', spkgrp, 'intens', intens, 'saveVar', true,...
    'force', true, 'extension', 'dat', 'recSystem', 'oe',...
    'protocol', 'io', 'anaflag', true, 'inspect', false, 'fsIn', fs,...
    'cf', 0);  

intens = (100);
fepsp = fEPSPfromWCP('basepath', basepath, 'sfiles', [],...
    'sufx', 'stp_3pulse', 'force', true, 'protocol', 'stp',...
    'intens', intens, 'inspect', false, 'fs', 30000);

fepsp = fEPSP_analysis('fepsp', fepsp, 'basepath', basepath,...
    'force', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spike sorting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in this pipeline waveforms are extracted from the whitened data and are
% used all the way until manual curation (including calculation of
% isolation distance). Afterwards, cell explorer resnips the waveform for
% the dat file and detrends them such that the cell class is determined
% from the raw waveforms.

% spike detection from temp_wh
[spktimes, ~] = spktimesWh('basepath', basepath, 'fs', fs, 'nchans', nchans,...
    'spkgrp', spkgrp, 'saveVar', true, 'saveWh', true,...
    'graphics', false, 'force', true);

% spike rate
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
% dur = -420;
% t = [];
% nsClip('dur', dur, 't', t, 'bkup', true, 'grp', [3 : 4]);

% clean clusters after sorting 
cleanCluByFet('basepath', pwd, 'manCur', true, 'grp', [1 : 4])

% cut spk from dat and realign
fixSpkAndRes('grp', 3, 'dt', 0, 'stdFactor', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cell explorer metrics
cell_metrics = ProcessCellMetrics('session', session,...
    'manualAdjustMonoSyn', false, 'summaryFigures', false,...
    'debugMode', true, 'transferFilesFromClusterpath', false,...
    'submitToDatabase', false, 'getWaveformsFromDat', true);
% cell_metrics = CellExplorer('basepath', basepath);

% cluster validation
load([basename, '.spikes.cellinfo.mat'])
mu = [];
spikes = cluVal('spikes', spikes, 'basepath', basepath, 'saveVar', true,...
    'saveFig', false, 'force', true, 'mu', mu, 'graphics', false,...
    'vis', 'on', 'spkgrp', spkgrp);

% pyr vs. int
cc = cellclass('basepath', basepath,...
    'waves', cat(1, spikes.rawWaveform{:})', 'saveVar', true,...
    'graphics', false, 'fs', fs);

% firing rate
binsize = 60;
winBL = [0 300 * 60];
winCalc = [0 Inf];
fr = firingRate(spikes.times, 'basepath', basepath, 'graphics', false, 'saveFig', false,...
    'binsize', binsize, 'saveVar', true, 'smet', 'MA', 'winBL',...
    winBL, 'winCalc', winCalc);

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
lfpInterval = [150 * 60, 180 * 60];
lfp = getLFP('basepath', basepath, 'ch', [spkgrp{:}], 'chavg', {},...
    'fs', 1250, 'interval', lfpInterval, 'extension', 'lfp',...
    'savevar', true, 'forceL', true, 'basename', '');

% remove 50 Hz from signal
emgOrig = filterLFP(emgOrig, 'fs', 1250, 'stopband', [49 51],...
    'dataOnly', true, 'saveVar', false, 'graphics', false);

acc = EMGfromACC('basepath', basepath, 'fname', [basename, '.lfp'],...
    'nchans', nchans, 'ch', nchans - 2 : nchans, 'saveVar', true, 'fsIn', 1250,...
    'graphics', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sleep states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prep signal
[EMG, EEG, sigInfo] = as_prepSig([basename, '.lfp'], acc.mag,...
    'eegCh', [4 : 7], 'emgCh', [], 'saveVar', true, 'emgNchans', [],...
    'eegNchans', nchans, 'inspectSig', false, 'forceLoad', true,...
    'eegFs', 1250, 'emgFs', 1250, 'emgCf', [10 600]);

% manually create labels
labelsmanfile = [basename, '.AccuSleep_labelsMan.mat'];
AccuSleep_viewer(EEG, EMG, 1250, 1, [], labelsmanfile)

% classify with a network
netfile = 'D:\Code\slutskycode\extracellular in vivo\lfp\SleepStates\AccuSleep\trainedNetworks\net_210709_142246.mat';
ss = as_wrapper(EEG, EMG, sigInfo, 'basepath', basepath, 'calfile', [],...
    'viaGui', false, 'forceCalibrate', false, 'inspectLabels', false,...
    'saveVar', true, 'forceAnalyze', true, 'fs', 1250, 'netfile', netfile,...
    'graphics', true);

% inspect separation after classifying / manual scoring
as_stateSeparation(EEG, EMG, labels)

% calc psd in states
[psdStates, faxis, emgRMS] = psd_states('eeg', EEG, 'emg', EMG,...
    'labels', labels, 'fs', 1250, 'graphics', true);

% get confusion matrix between two labels
[ss.netPrecision, ss.netRecall] = as_cm(labels1, labels2);

% show only x hours of data
x = 11;
tidx = [1 : x * 60 * 60 * 1250];
lidx = [1 : x * 60 * 60];
AccuSleep_viewer(EEG(tidx), EMG(tidx), 1250, 1, labels(lidx), [])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% handle dat  files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% cat dat
nchans = 2;
fs = 3051.7578125;
newpath = mousepath;
datFiles{1} = 'G:\lh86\lh86_210228_070000\lh86_210228_070000.emg.dat';
datFiles{2} = 'G:\lh86\lh86_210228_190000\lh86_210228_190000.emg.dat';
sigInfo = dir(datFiles{1});
nsamps = floor(sigInfo.bytes / class2bytes('int16') / nchans);
parts{1} = round([195360091 722703787] / 24414.06 * fs);
parts{2} = round([1 722703787] / 24414.06 * fs);

catDatMemmap('datFiles', datFiles, 'newpath', newpath, 'parts', parts,...
    'nchans', nchans, 'saveVar', true)


% preproc dat
clip = [0, 10 * 60 * 1250];
datInfo = preprocDat('basepath', pwd,...
    'fname', [basename, '.lfp'], 'mapch', 1 : 19,...
    'rmvch', [], 'nchans', 19, 'saveVar', false, 'clip', clip,...
    'chunksize', 5e6, 'precision', 'int16', 'bkup', true);


