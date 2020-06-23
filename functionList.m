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
% IO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = 'H:\Data\Dat\lh46\lh46_200226';
store = 'Raw1';
fs = 24414.06;
blocks = [1, 2];
chunksize = 300;
mapch = [1 : 16];
mapch = [1 : 2 : 7, 2 : 2 : 8, 9 : 2 : 15, 10 : 2 : 16];
rmvch = [];
clip = cell(1, 1);

% tank to dat
[datInfo] = tdt2dat('basepath', basepath, 'store', store, 'blocks',  blocks,...
    'chunksize', chunksize, 'mapch', mapch, 'rmvch', rmvch, 'clip', clip);

% open ephy to dat
exp = [10];
rec = cell(max(exp), 1);
datInfo = preprocOE('basepath', basepath, 'exp', exp, 'rec', rec,...
    'rmvch', rmvch, 'mapch', mapch, 'concat', true, 'nchans', 35);

% pre-process dat (remove channels, reorder, etc.)
datInfo = preprocDat('basepath', basepath, 'fname', '', 'mapch', mapch,...
    'rmvch', rmvch, 'nchans', 35, 'saveVar', true,...
    'chunksize', 5e5, 'precision', 'int16', 'bkup', true);

% session info (cell explorer)
session = sessionTemplate(pwd, 'showGUI', true);
basepath = session.general.basePath;
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;

% acceleration
newch = length(mapch) - length(rmvch);
chAcc = [newch : -1 : newch - 2];
newpath = fileparts(datInfo.newFile);
getACCfromDat('basepath', newpath, 'fname', '',...
    'nchans', newch, 'ch', chAcc, 'force', false, 'saveVar', true,...
    'graphics', false, 'fs', 1250);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LFP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fEPSP from open ephys
intens = [20 40 60 80 100];
fepsp = getfEPSPfromOE('basepath', basepath, 'fname', '', 'nchans', 31,...
    'tet', spkgrp, 'intens', intens, 'concat', false, 'saveVar', true,...
    'force', true);  

% fEPSP from TDT
lfp = getfEPSPfromTDT('basepath', basepath, 'andname', 'stability2',...
    'blocks', blocks, 'mapch', mapch, 'ch', ch, 'clip', clip,...
    'saveVar', true, 'fdur', 0.1, 'concat', true);

% load lfp
lfp = getLFP('basepath', basepath, 'ch', [spkpgrp{:}], 'chavg', {},...
    'fs', 1250, 'interval', [0 inf], 'extension', 'lfp',...
    'savevar', true, 'forceL', true, 'basename', '');

% anesthesia states (see also aneStates_wrp)
[bs, iis, ep] = aneStates('ch', 1, 'basepath', basepath,...
    'basename', basename, 'graphics', graphics,...
    'saveVar', saveVar, 'saveFig', saveFig, 'forceA', forceA,...
    'binsize', 30, 'smf', 7, 'thrMet', 1);
        
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

% EMG from LFP
emglfp = getEMGfromLFP(double(lfp.data(:, :)),...
    'emgFs', 10, 'saveVar', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load
spikes = getSpikes('basepath', basepath, 'saveMat', true,...
    'noPrompts', true, 'forceL', false);

% review clusters
mu = find(spikes.isi > 1);
mu = sort([mu', 2, 4, 19]);
spikes = cluVal(spikes, 'basepath', basepath, 'saveVar', false,...
    'saveFig', false, 'force', true, 'mu', mu, 'graphics', true,...
    'vis', 'on');

% separation of SU and MU
plotIsolation(basepath, spikes, false)

% CCG
% low res
binSize = 0.001; dur = 0.05; % low res
binSize = 0.0001; dur = 0.02; % high res
[ccg t] = CCG({spikes.times{:}}, [], 'duration', dur, 'binSize', binSize);

u = sort([20 27]);
plotCCG('ccg', ccg(:, u, u), 't', t, 'basepath', basepath,...
    'saveFig', false, 'c', {'k'}, 'u', spikes.UID(u));
    
% cell classification
basepath = 'E:\Data\Dat\Refaela\120520_Bsl';
spikes = getSpikes('basepath', basepath, 'saveMat', false,...
    'noPrompts', true, 'forceL', false);

waves = spikes.maxwv';
cc = cellclass('waves', waves,...
    'fs', spikes.samplingRate, 'man', false, 'spikes', spikes); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 6: firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
winBL = [info.lns(1) info.lns(3)] * spikes.samplingRate * 60;
fr = FR(x, 'basepath', basepath, 'graphics', false, 'saveFig', false,...
    'binsize', 60, 'saveVar', false, 'smet', 'MA');

% unite all units
x = {sort(vertcat(spikes.times{:}))};
[fr.strd, ~, fr.tstamps] = calcFR(x, 'binsize', 60,...
    'winCalc', [1 Inf], 'smet', 'none');

filename = bz_BasenameFromBasepath(basepath);
filename = [filename '.Raw1.Info.mat'];
load(filename);
info.labels = {'BL', '2', '3', '4'};
lns = cumsum(info.blockduration / 60);
lns = [1e-6 lns];
info.lns = lns(1);
save(filename, 'info');

info.lns = cumsum(datInfo.nsamps / 20000 / 60);

f = figure;
subplot(2, 1, 1)
plotFRtime('fr', fr, 'units', true, 'spktime', spikes.times,...
    'avg', false, 'lns', info.lns, 'lbs', info.labels,...
    'raster', false, 'saveFig', false);
title('')
xlabel('')
subplot(2, 1, 2)
plotFRtime('fr', fr, 'units', false, 'spktime', spikes.times,...
    'avg', true, 'lns', info.lns, 'lbs', info.labels,...
    'raster', false, 'saveFig', false);
xlabel('Time [m]')
title('')
figname = 'Firing Rate';
export_fig(figname, '-tif', '-transparent')
        
[nunits, nbins] = size(fr.strd);
tFR = ([1 : nbins] / (60 / fr.binsize) / 60);
p = plot(tFR, log10(fr.strd));
hold on
plot(tFR, mean(log10(fr.strd)), 'k', 'LineWidth', 3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 7: concatenate spikes from different sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parentdir = 'E:\Data\Others\Buzsaki';
basepath = parentdir;
structname = 'spikes';
spikes = catStruct(parentdir, structname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 8: get video projection from ToxTrack file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'TestProject';
vid = getVid(filename, 'basepath', basepath, 'graphics', true, 'saveFig', false, 'saveVar', false);


