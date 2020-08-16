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
basepath = 'D:\VMs\shared\lh49\lh49_200319';
store = 'Raw1';
fs = 24414.06;
blocks = [11, 13, 18];
chunksize = 300;
mapch = [1 : 16];
% mapch = [1 : 2 : 7, 2 : 2 : 8, 9 : 2 : 15, 10 : 2 : 16];
rmvch = [];
clip = cell(1, 1);

% tank to dat
[datInfo] = tdt2dat('basepath', basepath, 'store', store, 'blocks',  blocks,...
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
    'rmvch', rmvch, 'nchans', 35, 'saveVar', true,...
    'chunksize', 5e5, 'precision', 'int16', 'bkup', true);

% session info (cell explorer)
session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'force', true, 'saveVar', true);
basepath = session.general.basePath;
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;

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

% create LFP
bz_LFPfromDat(filepath, 'noPrompts', true)

% fEPSP from open ephys
intens = [100 50 150 200 250 300];
fepsp = getfEPSPfromOE('basepath', char(filepath), 'fname', '', 'nchans', 27,...
    'spkgrp', f{i}.spkgrp, 'intens', intens, 'concat', false, 'saveVar', true,...
    'force', true, 'vis', 'off');

intens = [100 150 200 250 300];
fepsp = getfEPSPfromOE('basepath', basepath, 'fname', '', 'nchans', nchans,...
    'spkgrp', spkgrp, 'intens', intens, 'concat', false, 'saveVar', true,...
    'force', true);  

% fEPSP from TDT
lfp = getfEPSPfromTDT('basepath', basepath, 'andname', 'stability2',...
    'blocks', blocks, 'mapch', mapch, 'ch', ch, 'clip', clip,...
    'saveVar', true, 'fdur', 0.1, 'concat', true);

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
binSize = 0.001; dur = 0.05; % low res
binSize = 0.0001; dur = 0.02; % high res
[ccg t] = CCG({spikes.times{:}}, [], 'duration', dur, 'binSize', binSize);
u = sort([20 27]);
plotCCG('ccg', ccg(:, u, u), 't', t, 'basepath', basepath,...
    'saveFig', false, 'c', {'k'}, 'u', spikes.UID(u));
    
% cell classification
waves = spikes.maxwv';
cc = cellclass('waves', waves,...
    'fs', spikes.samplingRate, 'man', false, 'spikes', spikes); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 6: firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
winBL = [info.lns(1) info.lns(3)] * spikes.samplingRate * 60;
fr = FR(x, 'basepath', basepath, 'graphics', false, 'saveFig', false,...
    'binsize', 60, 'saveVar', false, 'smet', 'MA');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis across sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fEPSP_sessions
fr_sessions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% behavior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'TestProject';
vid = getVid(filename, 'basepath', basepath, 'graphics', true,...
    'saveFig', false, 'saveVar', false);


