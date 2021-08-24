function ops = opsKS(varargin)

% kilosort configuration. based on configFile384.m,
% make_eMouseChannelMap.m, and KilosortConfiguration (kilosortWrapper).
% creates channel map file according to spike groups (designed for
% tetrodes).
% 
% INPUT:
%   basepath    string. path to recording folder {pwd}. if multiple dat
%               files exist only the first will be processed
%   procpath    string. path to where processing will occur.
%               should be a fast ssd
%   fs          numeric. sampling frequency [hz]{20000}
%   nchans      numeric. number of channels in dat file.
%   spkgrp      array where each cell is the electrodes for a spike group. 
%   trange      numeric. time range to sort [s]{[0 Inf]}
%
% DEPENDENCIES
%
% TO DO LIST:
% 
% 08 nov 20 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'procpath', '');
addOptional(p, 'fs', 20000, @isnumeric);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'spkgrp', {}, @iscell);
addOptional(p, 'trange', [0 Inf], @isnumeric);

parse(p, varargin{:})
basepath    = p.Results.basepath;
procpath    = p.Results.procpath;
fs          = p.Results.fs;
nchans      = p.Results.nchans;
spkgrp      = p.Results.spkgrp;
trange      = p.Results.trange;

if isempty(procpath)
    procpath = basepath;
end
kcoords = []; 
ycoords = []; 
xcoords = [];
xrep = [20 40 60 80];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channel map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% list of channel indices (including dead \ non-ephys channels)
chanMap = [spkgrp{:}];
chanMap0ind = chanMap - 1;
% the first thing Kilosort does is reorder the data with data = data(chanMap, :).
% Now we declare which channels are "connected" in this normal ordering,
% meaning not dead or used for non-ephys data
connected = true(nchans, 1);
% now we define the horizontal (xcoords) and vertical (ycoords) coordinates
% of these channels. For dead or nonephys channels the values won't matter.
% These are in um here, but the absolute scaling doesn't really matter in
% the algorithm. In addition, often multi-shank probes or tetrodes will be
% organized into groups of channels that cannot possibly share spikes with
% the rest of the probe. This helps the algorithm discard noisy templates
% shared across groups. In this case, we set kcoords to indicate which
% group the channel belongs to.
badch = setdiff([1 : nchans], chanMap);
for igrp = 1 : length(spkgrp)
    nchGrp = length(spkgrp{igrp});
    kcoords = [kcoords, ones(1, nchGrp) * igrp];
    xcoords = [xcoords, xrep(1 : nchGrp)];
    ycoords = [ycoords, ones(1, nchGrp) * igrp * 20];
end
% xcoords(badch) = NaN;
% ycoords(badch) = NaN;
% kcoords(badch) = NaN;
connected(badch) = false; % e.g. acceleration
% at this point in Kilosort we do data = data(connected, :), ycoords =
% ycoords(connected), xcoords = xcoords(connected) and kcoords =
% kcoords(connected) and no more channel map information is needed (in particular
% no "adjacency graphs" like in KlustaKwik).
% Now we can save our channel map and also a channel_shanks file for phy.
ops.chanMap = fullfile(basepath, 'chanMap.mat');
save(ops.chanMap, 'chanMap', 'chanMap0ind', 'connected',...
    'xcoords', 'ycoords', 'kcoords', 'fs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ks parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% arguments
ops.kcoords         = kcoords;
ops.fs              = fs;       % sampling rate
ops.NchanTOT        = nchans;   % total number of channels in your recording
ops.trange          = trange;   % time range to sort [s]

% kcoords is used to forcefully restrict templates to channels in the same
% channel group. criterionNoiseChannels can be set to allow a fraction of
% all templates to span more channel groups, so that they can capture
% shared noise across all channels. if this number is less than 1, it will
% be treated as a fraction of the total number of clusters. if this number
% is larger than 1, it will be treated as the "effective number" of channel
% groups at which to set the threshold. So if a template occupies more than
% this many channel groups, it will not be restricted to a single channel
% group {0.2}.
ops.criterionNoiseChannels = 0.2; 

% Thresholds on spike detection used during the optimization Th(1) or
% during the final pass Th(2). These thresholds are applied to the template
% projections, not to the voltage. Typically, Th(1) is high enough that the
% algorithm only picks up sortable units, while Th(2) is low enough that it
% can pick all of the spikes of these units. It doesn't matter if the final
% pass also collects noise: an additional per neuron threshold is set
% afterwards, and a splitting step ensures clusters with multiple units get
% split. {[10 4]}
ops.Th = [8 4];

% how important is the amplitude penalty {10}. The individual spike
% amplitudes are biased towards the mean of the cluster by this factor; 50
% is a lot, 0 is no bias.
ops.lam = 25;

% Threshold on the area under the curve (AUC) criterion for performing a
% split in the final step. If the AUC of the split is higher than this,
% that split is considered good. However, a good split only goes through
% if, additionally, the cross-correlogram of the split units does not
% contain a big dip at time 0. splitting a cluster at the end requires at
% least this much isolation for each sub-cluster (max = 1){0.9}.
ops.AUCsplit = 0.86;

% frequency for high pass filtering {150}
ops.fshigh = 300;

% frequency for low (band) pass filtering {none}
% ops.fslow = 300;

% minimum firing rate on a "good" channel (0 to skip){0.1}
ops.minfr_goodchannels = 0.05;

% minimum spike rate (Hz), if a cluster falls below this for too long it
% gets removed {1/50}.
ops.minFR = 1/100;

% number of samples to average over (annealed from first to second value)
% {20 400}
ops.momentum = [20 800];

% spatial constant in um for computing residual variance of spike
ops.sigmaMask = 30;

% spatial scale for datashift kernel (new in v2.5)
ops.sig = 20;

% type of data shifting (0 = none, 1 = rigid, 2 = nonrigid)
ops.nblocks = 5;

% threshold crossings for pre-clustering (in PCA projection space)
ops.ThPre = 8;

% number of time samples for the templates (has to be <=81 due to GPU
% shared memory). {61}. 
ops.nt0             = 61;       

% time sample where the negative peak should be aligned {2} 
ops.nt0min  = ceil(20 * ops.nt0 / 61); 

ops.spkTh           = -4;       % spike threshold in standard deviations {-6}
ops.reorder         = 1;        % whether to reorder batches for drift correction.
ops.nskip           = 25;       % how many batches to skip for determining spike PCs
ops.Nfilt           = 256;      % number of filters to use (2-4 times more than Nchan, should be a multiple of 32)
ops.nfilt_factor    = 8;        % max number of clusters per good channel (even temporary ones) {4}
ops.ntbuff          = 64;       % samples of symmetrical buffer for whitening and spike detection
ops.nSkipCov        = 25;       % skip n batches when computing whitening matrix 
ops.scaleproc       = 200;      % int16 scaling of whitened data
ops.nPCs            = 3;        % how many PCs to project the spikes into
ops.GPU             = 1;        % has to be 1, no CPU version yet, sorry
ops.useRAM          = 0;        % not yet available
ops.fig             = 1;        % plot figures

% nomber of time points per batch (try decreasing if out of memory). must
% be multiple of 32 + ntbuff.
ops.NT              = 128 * 1024 + ops.ntbuff;

% how many channels to whiten together (Inf for whole probe whitening)
ops.whiteningRange  = min([64 sum(connected)]); 

% perform common average referencing by median {1}.
ops.CAR = 1;

datfile             = dir(fullfile(basepath, '*.dat')); % find the binary file
ops.fbinary         = fullfile(basepath, datfile(1).name);
ops.fproc           = fullfile(procpath, 'temp_wh.dat'); % path to processing file

end

% EOF