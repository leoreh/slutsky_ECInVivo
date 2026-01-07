function rez = ks_wrapper()

% wrapper for running kilosort. loads configuration, runs the algorithm and
% arranges the output in phy or neurosuite format.
%
% INPUT:
%   basepath    string. path to recording folder {pwd}. if multiple dat
%               files exist only the first will be processed
%   trange      numeric. time range to sort [s]{[0 Inf]}
%
% DEPENDENCIES
%   npy-matlab
%   kilosort3
%   ks_ops
%   ks2ns
%
% TO DO LIST:
%   # currently, only the gui can be used in cases where one channel of a
%   tetrode is bad (done)
%
% 22 may 20 LH      updates:
% 16 jun 20             outFormat
% 31 oct 20             updated to KS2.5
% 08 nov 20             split opsKS
% 01 jun 22             ks3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = pwd;
session = CE_sessionTemplate(basepath, 'viaGUI', false,...
    'forceDef', true, 'forceL', false, 'saveVar', true);
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;

ops = ks_ops('basepath', basepath, 'fs', fs, 'nchans', nchans,...
    'spkgrp', spkgrp, 'trange', [0 4 * 60 * 60]);
if ~exist(ops.fbinary, 'file')
    error('%s not found', ops.fbinary)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run ks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ops params for skipping first two functions
bytes       = get_file_size(ops.fbinary); % size in bytes of raw binary
nTimepoints = floor(bytes/ops.NchanTOT/2); % number of total timepoints
ops.tstart  = ceil(ops.trange(1) * ops.fs); % starting timepoint for processing data segment
ops.tend    = min(nTimepoints, ceil(ops.trange(2) * ops.fs)); % ending timepoint
ops.sampsToRead = ops.tend-ops.tstart; % total number of samples to read
ops.Nbatch      = ceil(ops.sampsToRead / ops.NT); % number of data batches
ops.Nchan = numel(ops.kcoords); % total number of good channels that we will spike sort
ops.Nfilt = getOr(ops, 'nfilt_factor', 4) * ops.Nchan; % upper bound on the number of templates we can have
rez.ops             = ops;
rez.temp.Nbatch     = rez.ops.Nbatch;
load('chanMap.mat')
rez.yc = ycoords';
rez.xc = xcoords';
rez.ycoords = rez.yc;
rez.xcoords = rez.xc;
dmin = median(diff(unique(rez.yc)));
fprintf('vertical pitch size is %d \n', dmin)
% The min and max of the y and x ranges of the channels
ymin = min(rez.yc);
ymax = max(rez.yc);
xmin = min(rez.xc);
xmax = max(rez.xc);
rez.ops.dmin = dmin;
rez.ops.yup = ymin:dmin/2:ymax; % centers of the upsampled y positions
% dminx = median(diff(unique(rez.xc)));
yunq = unique(rez.yc);
mxc = zeros(numel(yunq), 1);
for j = 1:numel(yunq)
    xc = rez.xc(rez.yc==yunq(j));
    if numel(xc)>1
       mxc(j) = median(diff(sort(xc))); 
    end
end
dminx = max(5, median(mxc));
rez.ops.dminx = dminx;
nx = round((xmax-xmin) / (dminx/2)) + 1;
rez.ops.xup = linspace(xmin, xmax, nx); % centers of the upsampled x positions
disp(rez.ops.xup) 


rez                 = preprocessDataSub(ops);       % preprocess data to create temp_wh.dat
rez                 = datashift2(rez, 0);           % shift data doesn't work for us (tetrodes?)
[rez, st3, tF]      = extract_spikes(rez);
rez                 = template_learning(rez, tF, st3);
[rez, st3, tF]      = trackAndSort(rez);
rez                 = final_clustering(rez, tF, st3);
rez                 = find_merges(rez, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output for manual curation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save rez
save(fullfile(basepath, 'rez.mat'), 'rez', '-v7.3');

% output to phy
kspath = fullfile(basepath, 'kilosort3');
mkdir(kspath)
rezToPhy2(rez, kspath);

% output to ns
ks2ns(rez)

end

% EOF



