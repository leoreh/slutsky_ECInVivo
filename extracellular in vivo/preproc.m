% this is :a wrapper for preprocessing extracellular data.
% contains calls to various functions.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basepath to recording folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = 'F:\Refaela\RA\DataTanks\Bruce\290519_1D';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: file conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
store = 'Raw1';
blocks = [1];
chunksize = 60;
mapch = [1, 2, 5, 6, 3, 4, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
rmvch = [9, 10, 11, 12, 13, 14, 15, 16];
clip{1} = [60 * 0.1 * 60 Inf];

% tank to dat
[info] = tdt2dat(basepath, 'Raw1', blocks, chunksize, mapch, rmvch, clip);

% ddt to dat
filenames{3} = 'block2_ddt_edit12H';
ddt2dat(basepath, mapch, rmvch, 'filenames', filenames)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: load LFP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chans = [1:16];
chavg = {1 : 4; 5 : 7; 8 : 11; 12 : 15};
chavg = {};
lfp = getLFP('basepath', basepath, 'filename', 'ec013.527.eeg', 'chans', chans, 'chavg', chavg, 'fs', 1250,...
    'interval', [0 inf], 'savevar', true);

%%% ripples
ripples = findRipples(lfp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: load EMG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% option 1:
blocks = [2];
rmvch = [2:4];
emg = getEMG(basepath, 'EMG1', blocks, rmvch);

% option 2:
emglfp = getEMGfromLFP(double(lfp.data(:, chans)), 'emgFs', 10, 'saveVar', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: load spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spikes = getSpikes('basepath', basepath, 'saveMat', true, 'noPrompts', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: review clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spikes = cluVal(basepath, spikes);

% compare number of spikes and clusters from clustering to curation
numSpikes = getNumSpikes(basepath, spikes);

% plot separation of SU and MU
plotIsolation(basepath, spikes, false)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4: CCH temporal dynamics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 5: cell classification based on waveform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CellClass = cellClass(cat(1, spikes.rawWaveform{:})', 'fs', 20000); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 6: calculate mean firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fr = FR(spikes.times, 'basepath', basepath, 'graphics', true, 'saveFig', false, 'saveVar', false);
fr = FR(spikes.times, 'basepath', basepath, 'graphics', true, 'saveFig', false, 'saveVar', false,...
    'winCalc', [0 25200]);

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


