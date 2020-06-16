% preproc_wrapper
basepath = 'D:\tempData\lh46_200226a';
cd(basepath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open ephys
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = 'H:\Data\Dat\lh52\2020-06-09_09-05-32';
rmvch = [10, 12, 13, 16, 17, 21, 23, 24];
mapch = [25 26 27 28 30 1 2 29 3 : 14 31 0 15 16 17 : 24 32 33 34] + 1;
exp = [10];
rec = cell(max(exp), 1);
datInfo = preprocOE('basepath', basepath, 'exp', exp, 'rec', rec,...
    'rmvch', rmvch, 'mapch', mapch, 'concat', true, 'nchans', 35);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% session info (cell explorer foramt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
session = sessionTemplate(pwd, 'showGUI', true);

% session params
basepath = session.general.basePath;
nchans = session.extracellular.nChannels;
ngrp = session.extracellular.nSpikeGroups;
badch = nchans : -1 : nchans - 2;
% badch = [];
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intens = [4 6 8 10];
fepsp = getfEPSPfromOE('basepath', basepath, 'fname', '', 'nchans', nchans,...
    'spkgrp', spkgrp, 'intens', intens, 'concat', false, 'saveVar', true,...
    'force', true);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kilosort
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare channel map
kcoords = [];
ycoords = [];
xcoords = [];
xrep = [20 40 60 80];
for i = 1 : ngrp
    l = length(spkgrp{i});
    kcoords = [kcoords, ones(1, l) * i];
    xcoords = [xcoords, xrep(1 : l)];
    ycoords = [ycoords, ones(1, l) * i * 20];
end
xcoords(badch) = NaN;
ycoords(badch) = NaN;
kcoords(badch) = NaN;

rez = runKS('basepath', basepath, 'fs', fs, 'nchans', nchans,...
    'badch', badch, 'ngrp', ngrp, 'kcoords', kcoords, 'ycoords', ycoords, 'xcoords', xcoords,...
    'saveFinal', true, 'viaGui', false, 'cleanDir', false, 'trange', [0 Inf]);

% load('chanMap.mat')
% [~, rez.ops.basename] = fileparts(basepath);
% rez.ops.root = basepath;
% rez.ops.savepath = basepath;
% rez.ops.basepath = basepath;
% rez.xcoords = xcoords;
% rez.ycoords = ycoords;
% rez.connected = connected;
% rez.ops.ForceMaxRAMforDat = 10000000000;
% rez.ops.parfor = true;
% ConvertKilosort2Neurosuite_KSW(rez)
% Kilosort2Neurosui te(rez)
ks2ns(rez)

% Phy2Neurosuite(basepath,basepath,'neurosuite')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cell explorer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% session info
session = sessionTemplate(pwd, 'showGUI', true);
% calculate
cell_metrics = ProcessCellMetrics('session', session);
% gui
cell_metrics = CellExplorer('metrics',cell_metrics); 
% multiples sessions
basepaths = {'/Volumes/buzsakilab/Buzsakilabspace/Datasets/GirardeauG/Rat08/Rat08-20130708',...
    '/Volumes/buzsakilab/Buzsakilabspace/Datasets/GirardeauG/Rat08/Rat08-20130709'};

cell_metrics = LoadCellMetricsBatch('basepaths',basepaths);
cell_metrics = CellExplorer('metrics',cell_metrics);
% load a subset of units fullfilling multiple of criterium
% Get cells that are have a tag 'InverseSpike' or 'Good' and are assigned as 'Interneuron'
cell_metrics_idxs2 = LoadCellMetrics('cell_metrics',cell_metrics,'tags',...
    {'InverseSpike','Good'},'putativeCellType',{'Interneuron'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearKSdir(basepath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
binsize = 60;
fr = FR(spikes.times, 'basepath', basepath, 'graphics', false, 'saveFig', false,...
    'binsize', binsize, 'saveVar', false, 'smet', 'MA', 'winBL', [20 50 * 60]);

% unite all units
x = {sort(vertcat(spikes.times{:}))};
[fr.strd, ~, fr.tstamps] = calcFR(x, 'binsize', 60,...
    'winCalc', [1 Inf], 'smet', 'none');

% subplot(2, 1, 1)
stdshade(fr.norm, 0.3, 'k', fr.tstamps / 60, 3)
% subplot(2, 1, 2)
% plot(acc.tband / 60, acc.pband)

lbs = {'BL', '2', '3', '4'};
info.lns = lns(1);
hold on
addLns('lns', lns, 'lbs', lbs, 'ax', 'x')


