% this is a wrapper for preprocessing extracellular data.
% contains calls to various functions.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basepath to recording folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = 'E:\Data\Dat\lh37';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: file conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
store = 'Raw1';
fs = 24414.0125;
blocks = [9 : 14, 17 : 22, 25 : 32, 35 : 40];
chunksize = 60;
mapch = [1 : 16];
rmvch = [];
clip{28} = [0 1550];
clear clip
clip{1} = [];

% tank to dat
[info] = tdt2dat('basepath', basepath, 'store', store, 'blocks',  blocks,...
    'chunksize', chunksize, 'mapch', mapch, 'rmvch', rmvch, 'clip', clip);

% ddt to dat
filenames{3} = 'block2_ddt_edit12H';
ddt2dat(basepath, mapch, rmvch, 'filenames', filenames)

% dat to mat
[~, filename] = fileparts(basepath);
start = 0;      % s
duration = Inf; % s
nChannels = 16;
mat = bz_LoadBinary([filename '.dat'], 'frequency', fs, 'start', start,...
    'duration', duration, 'nChannels', nChannels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: load LFP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chans = [1:16];
chavg = {1 : 4; 5 : 7; 8 : 11; 12 : 15};
chavg = {};
lfp = getLFP('basepath', basepath, 'chans', chans, 'chavg', chavg, 'fs', 1250,...
    'interval', [0 inf], 'savevar', true, 'force', false);

%%% LFP events
events = findLFPevents('lfp', lfp, 'emgThr', 0, 'basepath', basepath,...
    'ch', [16], 'preset', 'bs');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: load EMG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% option 1:
blocks = [2];
rmvch = [2:4];
emg = getEMG(basepath, 'Stim', blocks, rmvch);

% option 2:
emglfp = getEMGfromLFP(double(lfp.data(:, chans)), 'emgFs', 10, 'saveVar', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4: load spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spikes = getSpikes('basepath', basepath, 'saveMat', true, 'noPrompts', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: review clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = [];
spikes = cluVal(spikes, 'basepath', basepath, 'saveVar', true,...
    'saveFig', true, 'force', true, 'mu', mu, 'graphics', true);

% compare number of spikes and clusters from clustering to curation
numSpikes = getNumSpikes(basepath, spikes);

% plot separation of SU and MU
plotIsolation(basepath, spikes, false)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4: CCH temporal dynamics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use buzcode CCG
% low res
binSize = 0.001; % [s]
dur = 0.05;
[ccg t] = CCG({spikes.times{:}}, [], 'duration', dur, 'binSize', binSize);
% high res
binSize = 0.0001; 
dur = 0.02;
[ccg t] = CCG({spikes.times{:}}, [], 'duration', dur, 'binSize', binSize);

for i = 1 : nunits
    nspikes(i) = length(spikes.times{i});
end

u = spikes.UID(nspikes > 6300);
u(1) = [];

u = sort([20 27]);
plotCCG('ccg', ccg(:, u, u), 't', t, 'basepath', basepath,...
    'saveFig', false, 'c', {'k'}, 'u', spikes.UID(u));

uu = datasample(u, 7, 'replace', false);
plotCCG('ccg', ccg(:, uu, uu), 't', t, 'basepath', basepath,...
    'saveFig', false, 'c', {'k'}, 'u', spikes.UID(uu));

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 5: cell classification based on waveform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CellClass = cellClass(cat(1, spikes.rawWaveform{:})', 'fs', spikes.samplingRate, 'man', true); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 6: calculate mean firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fr = FR(spikes.times, 'basepath', basepath, 'graphics', false, 'saveFig', false,...
    'binsize', 120, 'saveVar', true, 'smet', 'MA');

[~, filename] = fileparts(basepath);
filename = [filename '.Raw1.Info.mat'];
load(filename);
select = [1 : 6];
info.labels = {'Baseline', 'uPSEM 5 mg/kg', '', 'CNO 5 mg/kg'};
lns = cumsum(info.blockduration / 60 / 60);
lns = [1e-6, lns];
lns(end) = [];
lns = lns(select);
save(filename, 'info');

f = figure;
subplot(3, 1, 1)
plotFRtime('fr', fr, 'units', false, 'spktime', spikes.times,...
    'avg', false, 'lns', lns, 'lbs', info.labels(select),...
    'raster', true, 'saveFig', false);
subplot(3, 1, 2)
plotFRtime('fr', fr, 'units', true, 'spktime', spikes.times,...
    'avg', false, 'lns', lns, 'lbs', info.labels(select),...
    'raster', false, 'saveFig', false);
subplot(3, 1, 3)
plotFRtime('fr', fr, 'units', false, 'spktime', spikes.times,...
    'avg', true, 'lns', lns, 'lbs', info.labels(select),...
    'raster', false, 'saveFig', false);

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


