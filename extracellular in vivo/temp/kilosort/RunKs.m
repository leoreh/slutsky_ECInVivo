function rez = runKS(varargin)

% wrapper for running kilosort. based on configFile384.m,
% make_eMouseChannelMap.m, and KilosortConfiguration (kilosortWrapper).
% creates channel map file according to spike groups (designed for
% tetrodes). handles ks parameters arranged here according to those that
% should be played with. alas runs the algorithm and arranges the output in
% phy or neurosuite format.
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
%   saveFinal   logical. save final rez file {false}
%   viaGui      logical. run ks via gui (true) or script {false}
%   outFormat   string. output results to 'phy' or neurosuite {'ns'} format
%
% DEPENDENCIES
%   npy-matlab
%   kilosort 2
%   ks2ns
%
% TO DO LIST:
%   # currently, only the gui can be used in cases where one channel of a
%   tetrode is bad (done - 11 jun 20)
% 
% 22 may 20 LH      updates:
% 16 jun 20 LH      outFormat

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
addOptional(p, 'saveFinal', false, @islogical);
addOptional(p, 'viaGui', false, @islogical);
addOptional(p, 'cleanDir', false, @islogical);
addOptional(p, 'outFormat', 'ns', @ischar);

parse(p, varargin{:})
basepath    = p.Results.basepath;
procpath    = p.Results.procpath;
fs          = p.Results.fs;
nchans      = p.Results.nchans;
spkgrp      = p.Results.spkgrp;
trange      = p.Results.trange;
saveFinal   = p.Results.saveFinal;
viaGui      = p.Results.viaGui;
outFormat   = p.Results.outFormat;

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
chanMap = [1 : nchans];
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
badch = setdiff(chanMap, sort([spkgrp{:}]));
for i = 1 : length(spkgrp)
    l = length(spkgrp{i});
    kcoords = [kcoords, ones(1, l) * i];
    xcoords = [xcoords, xrep(1 : l)];
    ycoords = [ycoords, ones(1, l) * i * 20];
end
xcoords(badch) = NaN;
ycoords(badch) = NaN;
kcoords(badch) = NaN;
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

% Thresholds on spike detection used during the optimization Th(1) or
% during the final pass Th(2). These thresholds are applied to the template
% projections, not to the voltage. Typically, Th(1) is high enough that the
% algorithm only picks up sortable units, while Th(2) is low enough that it
% can pick all of the spikes of these units. It doesn't matter if the final
% pass also collects noise: an additional per neuron threshold is set
% afterwards, and a splitting step ensures clusters with multiple units get
% split. {[10 4]}
ops.Th = [10 2];

% how important is the amplitude penalty {10}. The individual spike
% amplitudes are biased towards the mean of the cluster by this factor; 50
% is a lot, 0 is no bias.
ops.lam = 12;

% Threshold on the area under the curve (AUC) criterion for performing a
% split in the final step. If the AUC of the split is higher than this,
% that split is considered good. However, a good split only goes through
% if, additionally, the cross-correlogram of the split units does not
% contain a big dip at time 0. splitting a cluster at the end requires at
% least this much isolation for each sub-cluster (max = 1){0.9}.
ops.AUCsplit = 0.88;

% frequency for high pass filtering {150}
ops.fshigh = 500;

% minimum firing rate on a "good" channel (0 to skip){0.1}
ops.minfr_goodchannels = 0;

% minimum spike rate (Hz), if a cluster falls below this for too long it
% gets removed {1/50}.
ops.minFR = 1/100;

% number of samples to average over (annealed from first to second value)
% {20 400}
ops.momentum = [20 800];

% spatial constant in um for computing residual variance of spike
ops.sigmaMask = 30;

% threshold crossings for pre-clustering (in PCA projection space)
ops.ThPre = 8;

ops.nt0             = 61;       % window width in samples. 2 ms for 20 kH {61}
ops.spkTh           = -4;       % spike threshold in standard deviations {-6}
ops.reorder         = 1;        % whether to reorder batches for drift correction.
ops.nskip           = 25;       % how many batches to skip for determining spike PCs
ops.Nfilt           = 256;      % number of filters to use (2-4 times more than Nchan, should be a multiple of 32)
ops.nfilt_factor    = 8;        % max number of clusters per good channel (even temporary ones)
ops.ntbuff          = 64;       % samples of symmetrical buffer for whitening and spike detection
ops.nSkipCov        = 25;       % compute whitening matrix from every N-th batch
ops.scaleproc       = 200;      % int16 scaling of whitened data
ops.nPCs            = 3;        % how many PCs to project the spikes into
ops.GPU             = 1;        % has to be 1, no CPU version yet, sorry
ops.useRAM          = 0;        % not yet available

% batch size (try decreasing if out of memory). must be multiple of 32 + ntbuff.
ops.NT              = 128 * 1024 + ops.ntbuff;

% how many channels to whiten together (Inf for whole probe whitening)
ops.whiteningRange  = min([64 sum(connected)]); 

datfile             = dir(fullfile(basepath, '*.dat')); % find the binary file
ops.fbinary         = fullfile(basepath, datfile(1).name);
ops.fproc           = fullfile(procpath, 'temp_wh.dat'); % path to processing file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run ks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Looking for data inside %s \n', basepath)

if viaGui
    kilosort
else
    % preprocess data to create temp_wh.dat
    rez = preprocessDataSub(ops);
    
    % time-reordering as a function of drift
    rez = clusterSingleBatches(rez);
    
    % intermediate saving because the rest can be resumed after loading rez
    save(fullfile(basepath, 'rez.mat'), 'rez', '-v7.3');
    
    % main tracking and template matching algorithm
    rez = learnAndSolve8b(rez);
    
    % final merges
    rez = find_merges(rez, 1);
    
    % final splits by SVD
    rez = splitAllClusters(rez, 1);
    
    % final splits by amplitudes
    rez = splitAllClusters(rez, 0);
    
    % decide on cutoff
    rez = set_cutoff(rez);
    
    fprintf('found %d good units \n', sum(rez.good > 0))
    
    % save final results
    if saveFinal
        rez = rmfield(rez, 'st2');
        rez = rmfield(rez, 'cProjPC');
        rez = rmfield(rez, 'cProj');
        fprintf('\nsaving final rez file\n')
        save(fullfile(basepath, 'rez.mat'), 'rez', '-v7.3');
    end
    
    % write for manual curation
    switch outFormat
        case 'phy'
            fprintf('\nSaving results to Phy \n')
            % note this function delets all .npy in path
            rezToPhy(rez, basepath); 
            % save channel_shanks file
            writeNPY(kcoords(~isnan(kcoords)), fullfile(basepath, 'channel_shanks.npy'));
        case 'ns'
            fprintf('\nSaving results to Neurosuite format \n')
            ks2ns(rez)
    end
    
    % clean temp_wh file
    delete(fullfile(basepath, 'temp_wh.dat'))
end



