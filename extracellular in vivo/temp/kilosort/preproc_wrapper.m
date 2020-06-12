% preproc_wrapper


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open ephys
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = 'H:\Data\Dat\lh52\2020-06-11_08-58-30';
cd(basepath)
rmvch = [10, 12, 13, 16, 17, 21, 23, 24];
mapch = [25 26 27 28 30 1 2 29 3 : 14 31 0 15 16 17 : 24 32 33 34] + 1;
exp = [3, 4];
rec = cell(max(exp), 1);
datInfo = preprocOE('basepath', basepath, 'exp', exp, 'rec', rec,...
    'rmvch', rmvch, 'mapch', mapch, 'concat', true, 'nchans', 35);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% session info (cell explorer foramt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
session = sessionTemplate(pwd, 'showGUI', true);

% session params
basepath = fileparts(datInfo.newFile);
nchans = session.extracellular.nChannels;
ngrp = session.extracellular.nSpikeGroups;
badch = nchans : -1 : nchans - 2;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intens = [4 6 8 10 12 15 20];
fepsp = getfEPSPfromOE('basepath', basepath, 'fname', '', 'nchans', nchans,...
    'spkgrp', spkgrp, 'intens', intens, 'concat', false, 'saveVar', true,...
    'force', true);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kilosort
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare channel map
kc = [];
yc = [];
xc = [];
xrep = [20 40 60 80];
for i = 1 : ngrp
    l = length(spkgrp{i});
    kc = [kc, ones(1, l) * i];
    xc = [xc, xrep(1 : l)];
    yc = [yc, ones(1, l) * i * 20];
end
xc(badch) = NaN;
yc(badch) = NaN;
kc(badch) = NaN;

rez = runKS('basepath', basepath, 'fs', fs, 'nchans', nchans,...
    'badch', badch, 'ngrp', ngrp, 'kc', kc, 'yc', yc, 'xc', xc,...
    'saveFinal', false, 'viaGui', false, 'cleanDir', false, 'trange', [0 Inf]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cell explorer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
