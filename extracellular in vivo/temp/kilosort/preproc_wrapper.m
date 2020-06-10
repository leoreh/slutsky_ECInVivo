% preproc_wrapper

basepath = 'F:\Data\Dat\lh52\lh52_200609\100642';
cd(basepath)

fs = 20000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kcoords = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 5 5 6 6 7 7 7 8 8];
kcoords_full = sort(repmat(1 : ngrp, 1, 4));
nchans = 27;
ngrp = 8;
badch = nchans : -1 : nchans - 2;

% list of channel indices (including dead \ non-ephys channels)
chanMap = [1 : nchans];
chanMap0ind = chanMap - 1;
% the first thing Kilosort does is reorder the data with data = data(chanMap, :).
% Now we declare which channels are "connected" in this normal ordering,
% meaning not dead or used for non-ephys data
connected = true(nchans, 1);
connected(badch) = false; % e.g. acceleration
% now we define the horizontal (x) and vertical (y) coordinates of these
% channels. For dead or nonephys channels the values won't matter. Again
% I will take this information from the specifications of the probe. These
% are in um here, but the absolute scaling doesn't really matter in the
% algorithm.
xcoords = repmat([20 40 60 80], 1, ngrp);
xcoords(badch) = NaN;
ycoords = sort(repmat([20 : 20 : ngrp * 20], 1, 4));
ycoords(badch) = NaN;
% Often, multi-shank probes or tetrodes will be organized into groups of
% channels that cannot possibly share spikes with the rest of the probe. This helps
% the algorithm discard noisy templates shared across groups. In
% this case, we set kcoords to indicate which group the channel belongs to.
kcoords = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 7 7 7 7];
kcoords(29 : 31) = NaN;
% at this point in Kilosort we do data = data(connected, :), ycoords =
% ycoords(connected), xcoords = xcoords(connected) and kcoords =
% kcoords(connected) and no more channel map information is needed (in particular
% no "adjacency graphs" like in KlustaKwik).
% Now we can save our channel map and also a channel_shanks file for phy.

rez = runKS('basepath', basepath, 'fs', fs, 'nchans', nchans,...
    'badch', badch, 'ngrp', ngrp, 'saveFinal', false, 'viaGui', false,...
    'cleanDir', false, 'trange', [0 Inf]);

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
