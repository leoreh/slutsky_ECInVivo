% preproc_wrapper
basepath = 'D:\VMs\shared\lh52_200625\170630';
cd(basepath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open ephys
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = 'E:\Leore\lh52\2020-06-21_09-07-33';
rmvch = [10, 12, 13, 16, 17, 21, 23, 24];
mapch = [25 26 27 28 30 1 2 29 3 : 14 31 0 15 16 17 : 24 32 33 34] + 1;
exp = [6];
rec = cell(max(exp), 1);
% rec{3} = 3;
intens = 
datInfo = preprocOE('basepath', basepath, 'exp', exp, 'rec', rec,...
    'rmvch', rmvch, 'mapch', mapch, 'concat', true, 'nchans', 35,...
    'intens', intens);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% session info (cell explorer foramt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
session = sessionTemplate(pwd, 'showGUI', true);
basepath = session.general.basePath;
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intens = [100 150 200 250 300];
fepsp = getfEPSPfromOE('basepath', basepath, 'fname', '', 'nchans', nchans,...
    'spkgrp', spkgrp, 'intens', intens, 'concat', false, 'saveVar', true,...
    'force', true);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kilosort
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rez = runKS('basepath', basepath, 'fs', fs, 'nchans', nchans,...
    'spkgrp', spkgrp, 'saveFinal', true, 'viaGui', false,...
    'cleanDir', false, 'trange', [0 Inf], 'outFormat', 'ns');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fix manual curation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fixSpkAndRes('grp', 6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cell explorer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load waveforms from spk files
% spikes = loadSpikes('session', session, 'useNeurosuiteWaveforms', true,...
%     'forceReload', true);
cell_metrics = ProcessCellMetrics('session', session,...
    'manualAdjustMonoSyn', false, 'summaryFigures', false,...
    'debugMode', true, 'transferFilesFromClusterpath', false,...
    'submitToDatabase', false);
cell_metrics = CellExplorer('metrics',cell_metrics); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearKSdir(basepath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% su / mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isi = cell_metrics.refractoryPeriodViolation; % percent isi < 2 ms
mu = find(isi < 10);
spikes = cluVal(spikes, 'basepath', basepath, 'saveVar', false,...
    'saveFig', false, 'force', true, 'mu', mu, 'graphics', false,...
    'vis', 'on', 'spkgrp', spkgrp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sleep states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badch = setdiff([session.extracellular.electrodeGroups.channels{:}],...
    [session.extracellular.spikeGroups.channels{:}]);
SleepScoreMaster(basepath, 'rejectChannels', badch)
TheStateEditor

ce_LFPfromDat(session)
bz_LFPfromDat(filepath)

EMGFromLFP = bz_EMGFromLFP(basepath, 'overwrite', false,...
    'rejectChannels', [25 : 27], 'noPrompts', false,...
    'saveMat', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
binsize = 60;
fr = FR(spikes.times, 'basepath', basepath, 'graphics', false, 'saveFig', false,...
    'binsize', binsize, 'saveVar', false, 'smet', 'MA', 'winBL', [20 50 * 60]);

% get pyr
pyr = strcmp(cell_metrics.putativeCellType, 'Pyramidal Cell');
int = strcmp(cell_metrics.putativeCellType, 'Narrow Interneuron');

figure
subplot(2, 1, 1)
stdshade(fr.strd(pyr, :), 0.3, 'k', fr.tstamps / 60, 3)
subplot(2, 1, 2)
stdshade(fr.strd(int, :), 0.3, 'k', fr.tstamps / 60, 3)



stdshade(fr.norm, 0.3, 'k', fr.tstamps / 60, 3)
% subplot(2, 1, 2)
% plot(acc.tband / 60, acc.pband)

lbs = {'BL', '2', '3', '4'};
info.lns = lns(1);
hold on
addLns('lns', lns, 'lbs', lbs, 'ax', 'x')


