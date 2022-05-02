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
    'forceDef', true, 'forceL', true, 'saveVar', true);      
basepath = session.general.basePath;
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;
[~, basename] = fileparts(basepath);

% add timebins to datInfo
[timebins, timepnt] = metaInfo_timebins('reqPnt', 5.5 * 60 * 60);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocessing of raw files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = 'I:\Data\lh107\2022-05-01_10-26-41';
mapch = 1 : 21;
rmvch = [18];
rmvch = [5, 17];
mapch = [1 : 6];
rmvch = [1, 3, 5];
mapch = [1 : 8];
rmvch = [2, 4, 6, 7, 8];

% tank to dat
store = 'Raw1';
blocks = [2 : 8];
chunksize = 300;
clip = cell(1, 1);
datInfo = tdt2dat('basepath', basepath, 'store', store, 'blocks',  blocks,...
    'chunksize', chunksize, 'mapch', mapch, 'rmvch', rmvch, 'clip', clip);

% open ephys to dat
exp = [1 : 3];
rec = cell(max(exp), 1);
datInfo = preprocOE('basepath', basepath, 'exp', exp, 'rec', rec,...
    'rmvch', rmvch, 'mapch', mapch,...
    'nchans', length(mapch), 'fsIn', 20000);

% digital input from OE
orig_paths = cellfun(@(x) fileparts(x), orig_files, 'uni', false);
orig_paths = cellfun(@(x) fileparts(x), orig_paths, 'uni', false);
orig_paths = cellfun(@(x) fileparts(x), orig_paths, 'uni', false);
exPathNew = pwd;
getDinOE('basepath', orig_paths, 'newpath', exPathNew,...
    'concat', true, 'saveVar', true);

% concatenate timestamps.npy and make sure dat files are not zero padded
cat_OE_tstamps('orig_paths', orig_paths, 'new_path', exPathNew,...
    'nchans', length(mapch));

% pre-process dat (remove channels, reorder, etc.)
clear orig_files
orig_files{1} = 'I:\Data\lh107\lh107_220430_105048\lh107_220430_105048.dat';
% orig_files{2} = 'L:\Data\lh101\lh101_220417_140957\lh101_220417_140957.dat';

clip{1} = [seconds(minutes(839)), seconds(minutes(882))]  * 20000;
% clip = [];

datInfo = preprocDat('orig_files', orig_files, 'mapch', mapch,...
    'rmvch', [rmvch], 'nchans', length(mapch), 'saveVar', true,...
    'chunksize', 1e7, 'precision', 'int16', 'clip', clip);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LFP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create lfp file
LFPfromDat('basepath', basepath, 'cf', 450, 'chunksize', 5e6,...
    'nchans', nchans, 'fsOut', 1250, 'fsIn', fs)    

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
    'graphics', false, 'force', true);

% get ripples
ripp = getRipples('basepath', basepath, 'rippCh', [9],...
    'emg', acc.mag, 'recWin', [0, Inf], 'saveVar', true,...
    'spkFlag', true, 'graphics', true, 'saveVar', true);

% remove pli
[EMG, tsaSig, ~] = tsa_filter('sig', EMG, 'fs', fs, 'tw', true,...
    'ma', true, 'graphics', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sleep states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load sSig or just the sig info
sigInfo = load([basename, '.sleep_sig.mat'], 'info');
sSig = load([basename, '.sleep_sig.mat']);

% call for acceleration
sSig = as_prepSig([basename, '.lfp'], acc.mag,...
    'eegCh', [1 : 4], 'emgCh', [], 'saveVar', true, 'emgNchans', [],...
    'eegNchans', nchans, 'inspectSig', true, 'forceLoad', true,...
    'eegFs', 1250, 'emgFs', 1250, 'eegCf', [], 'emgCf', [10 450], 'fs', 1250);

% call for lfp
sSig = as_prepSig([basename, '.lfp'], [],...
    'eegCh', [2], 'emgCh', [3], 'saveVar', true, 'emgNchans', [],...
    'eegNchans', nchans, 'inspectSig', true, 'forceLoad', true,...
    'eegFs', 1250, 'emgFs', 1250, 'eegCf', [], 'emgCf', [10 450], 'fs', 1250);

% manually create labels
labelsmanfile = [basename, '.sleep_labelsMan.mat'];
AccuSleep_viewer(sSig, [], labelsmanfile)

% classify with a network
netfile = 'D:\Code\slutsky_ECInVivo\lfp\SleepStates\AccuSleep\trainedNetworks\net_220107_233307.mat';
netfile = [];
% calData = ss.info.calibrationData;
calData = [];
ss = as_classify(sSig, 'basepath', pwd, 'inspectLabels', false,...
    'saveVar', false, 'forceA', true, 'netfile', netfile,...
    'graphics', true, 'calData', calData);

% inspect separation after classifying / manual scoring
as_stateSeparation(sSig, ss)

% get confusion matrix between two labels
[ss.netPrecision, ss.netRecall] = as_cm(labels1, labels2);

% plot state duration in timebins
[totDur, epLen] = as_plotZT('nwin', 8, 'sstates', [1, 2, 3, 4, 5],...
    'ss', ss, 'timebins', session.general.timebins);

% calc psd in states
psdBins = psd_states_timebins('basepath', pwd,...
    'chEeg', [], 'forceA', true, 'graphics', true,...
    'timepoints', timepoints, 'nbins', 8, 'saveVar', false);

% calc spec
spec = calc_spec('sig', [], 'fs', 1250, 'graphics', true, 'saveVar', true,...
    'padfft', -1, 'winstep', 1, 'logfreq', false, 'ftarget', [],...
    'ch', [{1 : 4}, {18}], 'force', true);

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
    'graphics', false, 'force', true, 'winWh', [0 Inf]);

% spike rate per tetrode. note that using firingRate requires
% special care becasue spktimes is given in samples and not seconds
for igrp = 1 : length(spkgrp)
    spktimes{igrp} = spktimes{igrp} / fs;
end
sr = firingRate(spktimes, 'basepath', basepath,...
    'graphics', true, 'binsize', 60, 'saveVar', 'sr', 'smet', 'none',...
    'winBL', [0 Inf]);

% create ns files 
dur = [];
t = [];
spktimes2ns('basepath', basepath, 'fs', fs,...
    'nchans', nchans, 'spkgrp', spkgrp, 'mkClu', true,...
    'dur', dur, 't', t, 'grps', [1 : length(spkgrp)],...
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
excludeMetrics = [{'deepSuperficial'},...
    {'theta_metrics'}, {'spatial_metrics'}, {'event_metrics'},...
    {'manipulation_metrics'}, {'state_metrics'}, {'psth_metrics'}];
cell_metrics = ProcessCellMetrics('session', session,...
    'manualAdjustMonoSyn', false, 'summaryFigures', false,...
    'debugMode', false, 'transferFilesFromClusterpath', false,...
    'submitToDatabase', false, 'getWaveformsFromDat', true,...
    'forceReloadSpikes', false, 'showFigures', false,...
    'excludeMetrics', excludeMetrics);
% cell_metrics = CellExplorer('basepath', pwd);

% firing rate
load([basename, '.spikes.cellinfo.mat'])
if isfield(session.general, 'timepnt')
    timepnt = session.general.timepnt;
else
    timepnt = Inf;
end
winBL = [0 timepnt];
fr = firingRate(spikes.times, 'basepath', basepath,...
    'graphics', true, 'binsize', 60, 'saveVar', true,...
    'smet', 'none', 'winBL', winBL, 'winCalc', [0, Inf]);

% select specific units
units = selectUnits('basepath', pwd, 'grp', [], 'saveVar', true,...
    'forceA', true, 'frBoundries', [0 Inf; 0 Inf],...
    'spikes', []);  

% plot fr vs. time
unitsClean = units.idx & units.gini' & units.mfrBL';
plot_FRtime_session('basepath', pwd,...
    'muFlag', false, 'saveFig', false,...
    'dataType', 'strd', 'units', units.idx)

% number of units per spike group
plot_nunits_session('basepath', basepath, 'frBoundries', [])

% cluster validation
spikes = cluVal('spikes', spikes, 'basepath', basepath, 'saveVar', true,...
    'saveFig', false, 'force', true, 'mu', [], 'graphics', false,...
    'vis', 'on', 'spkgrp', spkgrp);

% spike timing metrics
st = spktimes_metrics('winCalc', [], 'forceA', true);

% burstiness (mea)
brst = spktimes_meaBrst(spikes.times, 'binsize', Inf, 'isiThr', 0.02,...
    'minSpks', 2, 'saveVar', true, 'force', true);

% spike waveform metrics
swv = spkwv_metrics('basepath', basepath, 'fs', fs, 'forceA', true);

% mfr by states in time bins
timebins = session.general.timebins;
fr_timebins('basepath', pwd, 'forceA', true, 'graphics', false,...
    'timebins', timebins, 'saveVar', true, 'sstates', [1, 4, 5]);

% mono synaptic interactions (spike transmission gain)
monosyn = monoSyn_wrapper('spktimes', spikes.times, 'basepath', pwd,...
    'winCalc', [0, Inf], 'saveVar', true, 'graphics', true,...
    'forceA', true, 'fs', fs, 'saveFig', false,...
    'wv', swv.wv, 'wv_std', swv.wv_std);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mea
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname = '2021-05-16T17-34-298hBsln_20uMK_Bac_spikes SORTED';
mea = mea_orgNex('fname', fname, 'basepath', pwd, 'forceL', true);

mea_analyze('basepath', pwd,...
    'winBL', [0, 120 * 60], 'graphics', false, 'forceA', false)

% mono synaptic interactions (spike transmission gain)
monosyn = monoSyn_wrapper('spktimes', mea.spktimes, 'basepath', pwd,...
    'winCalc', [0, Inf], 'saveVar', true, 'graphics', true,...
    'forceA', true, 'fs', mea.info.fs, 'saveFig', false,...
    'wv', mea.wv, 'wv_std', mea.wv_std);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% utilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% concatenate var from different sessions
mname = 'lh104';
basepaths = {basepath};
[expData, xData] = sessions_catVarTime('mname', mname,...
    'dataPreset', {'emg', 'fr'}, 'graphics', true,...
    'basepaths', basepaths, 'xTicksBinsize', 1, 'markRecTrans', true);

% snip segments (e.g. spikes) from binary 
[spkwv, ~] = snipFromBinary('stamps', spktimes, 'fname', datname,...
    'win', win, 'nchans', nchans, 'ch', ch, 'align_peak', 'min',...
    'precision', 'int16', 'rmv_trend', 6, 'saveVar', false,...
    'l2norm', false);
        
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

% loads a list of variabels from multiple sessions and organizes them in a
% struct
varArray = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);

% set matlab graphics to custom or factory defaults
setMatlabGraphics(false)
