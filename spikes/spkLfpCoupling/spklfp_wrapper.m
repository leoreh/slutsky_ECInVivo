
% spklfp_wrapper

% investigate spk - lfp coupling.

% spklfp_singleband calculates for a specifiec band (e.g. theta) the phase coupling and rate
% map as a function of phase and lfp power. spklfp_wideband calculates the
% correlation between spike counts in 5 ms bins and the lfp magnitude and
% the correlation between population synchrony (in 5 ms bins) and lfp
% magnitude.

% in addition to working on multiple / single frequency, single band
% filters the data by fir whereas wide band filters each frequency band
% with by wavelet convolution. 

% based on bz_GenSpikeLFPCoupling, bz_PhaseModulation, and
% bz_PowerPhaseRatemap. alternatives include:
% (1) coherencypt from chronux which calculates the cross
% multitaper spectroms of the (point-process) spktimes and (contineous)
% lfp. see http://www-users.med.cornell.edu/~jdvicto/pdfs/pubo08.pdf 
% (2) the sta method of Vinck 2012 which is implemented in fieldtrip. see
% https://www.fieldtriptoolbox.org/tutorial/spikefield/ however, i think
% both of these methods are mainly suitable for tasks with repeated trials
% i think

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = pwd;
[~, basename] = fileparts(basepath);
cd(basepath)

% files
spksfile = fullfile(basepath, [basename, '.spikes.cellinfo.mat']);
unitsfile = fullfile(basepath, [basename, '.units.mat']);
lfpfile = fullfile(basepath, [basename, '.lfp']);
sessionfile = fullfile(basepath, [basename, '.session.mat']);
spklfpfile = fullfile(basepath, [basename, '.spklfp.mat']);

% load session info
if ~exist(sessionfile, 'file')
    session = CE_sessionTemplate(pwd, 'viaGUI', false,...
        'forceDef', false, 'forceL', false, 'saveVar', false);
else
    load(sessionfile)
end

% params from session info
nchans = session.extracellular.nChannels;
fs = session.extracellular.srLfp;

% load spikes
if exist(spksfile, 'file')
    load(spksfile, 'spikes')
end
nunits = length(spikes.times);

% load lfp
winCalc = [0 3 * 60 * 60];
recDur = diff(winCalc);
ch = 9 : 11;
sig = double(bz_LoadBinary(lfpfile, 'duration', diff(winCalc),...
    'frequency', fs, 'nchannels', nchans, 'start', winCalc(1),...
    'channels', ch, 'downsample', 1));
sig = mean(sig, 2);

% filter lfp


% clip spike times
spktimes = cellfun(@(x) x(InIntervals(x, winCalc)),...
    spikes.times, 'uni', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single band analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frange = [1 12];
[spklfp] = spklfp_singleband('basepath', basepath, 'fs', fs,...
    'sig', sig, 'spktimes', spktimes, 'frange', frange,...
    'jitterSig', false, 'srtUnits', true, 'graphics', true,...
    'saveVar', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wide band analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frange = [1 100];
nfreq = 20;
[spklfp] = spklfp_wideband('basepath', basepath, 'fs', fs,...
    'sig', sig, 'spktimes', spktimes, 'frange', frange,...
    'nfreq', nfreq, 'srtUnits', true, 'graphics', true,...
    'saveVar', true);

