function spec = calc_spec(varargin)

% calculates the multitaper spectrogram via chronux with refinements from
% accusleep. data can be loaded from .lfp or provided as input.
% multichannel support includes averaging channels (e.g. electrodes from
% the same tetrode) or avergaing spectrograms across channels and/or
% calculating the spectrogram separately for groups of channels. these are
% determined by the input 'ch'
%
% INPUT
%   basepath    char. fullpath to recording folder {pwd}
%   sig         lfp signal. if matrix assumes columns are samples and rows
%               are channels. if empty will load from (basename).lfp
%               according to ch
%   ch          cell of vecs depicting groups of channels (indices to the
%               rows of sig) whose spectrogram should be averaged. this is
%               done in mtspectrumc.m. if sig is loaded from binary then
%               averaging will be done on the lfp traces and not the
%               spctrogram to save computation time. this was done for
%               tetrodes (i.e. ch = spkgrp)
%   fs          sampling frequency {1250}.
%   ftarget     numeric. target frequency range and resolution. this can 
%               be used to control the frequency axis of the spectrogram
%               which depends on (1) the degree of zero padding chunks of
%               the signal for the fft ('pad') and (2) the time
%               resolution (frequency binsize = 1 / window). {[]}. if empty
%               than the freuqency range will be [0 120] and the resolution
%               will be determined by window and pad.
%   padfft      numeric. zero pad chunks of the data to the next pow2
%               when applying fft. 0 means yes -1 means no. see mtspecgramc.m
%   winstep     numeric. determines the time resolution of the spectrogram.
%               for accusleep should be equal to epoch length. {1} [sec]
%   logfreq     logical. ploy y axis (freq) on logscale {false}
%   graphics    logical. plot figure {false}
%   saveVar     logical. organize and save struct {true}
% 
% OUTPUT
%   spec        struct with fields:
%   s           spectrogram time x frequency x groups of channels.
%   
%
% CALLS
%   mtspecgramc
% 
% TO DO LIST
%       # find a way to set request freqs in a log scale
%       # normalize spectrogram
%       # separate graphics to stand alone (done)
%       # calc power in bands (e.g. delta / theta)
%
% 13 jan 20 LH      updates:
% 20 feb 20         normalize to broadband
% 25 feb 22         adapted mtspecgramc
% 18 apr 22         multichannel support

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar)
addParameter(p, 'sig', [], @isnumeric)
addParameter(p, 'ch', {}, @iscell)
addParameter(p, 'fs', 1250, @isnumeric)
addParameter(p, 'ftarget', [], @isnumeric)
addParameter(p, 'padfft', 0, @isnumeric)
addParameter(p, 'winstep', 1, @isnumeric)
addParameter(p, 'logfreq', false, @islogical)
addParameter(p, 'graphics', false, @islogical)
addParameter(p, 'saveVar', true, @islogical)

parse(p, varargin{:})
basepath        = p.Results.basepath;
sig             = p.Results.sig;
ch              = p.Results.ch;
fs              = p.Results.fs;
ftarget         = p.Results.ftarget;
padfft          = p.Results.padfft;
winstep         = p.Results.winstep;
logfreq         = p.Results.logfreq;
graphics        = p.Results.graphics;
saveVar         = p.Results.saveVar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% files
[~, basename] = fileparts(basepath);
specfile = fullfile(basepath, [basename, '.spec.mat']);
sessionfile = fullfile(basepath, [basename, '.session.mat']);

% session params
if exist(sessionfile, 'file')
    load(sessionfile, 'session')
    nchans = session.extracellular.nChannels;
    spkgrp = session.extracellular.spikeGroups.channels;
end

% prep frequencies
if isempty(ftarget)
    frange = [1 100];
else
    frange = [ftarget(1), ftarget(end)];
    if frange(2) > fs / 2
        error('requested max frequency greater than nyquist')
    end
end

% mtspecgramc params
window = max([5, winstep]);
mtspec_params.pad = padfft;
mtspec_params.Fs = fs;
mtspec_params.fpass = frange;
mtspec_params.tapers = [3 5];
mtspec_params.trialave = 1;

% organize channel groups
spec.info.ch = ch;
if isempty(ch)
    ch = {1};
end
ngrp = length(ch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prep sig
if isempty(sig)
    for igrp = 1 : ngrp
        sig(:, igrp) = mean(double(bz_LoadBinary([basename, '.lfp'], 'duration', Inf,...
            'frequency', 1250, 'nchannels', nchans, 'start', 0,...
            'channels', ch{igrp}, 'downsample', 1)), 2);
        ch{igrp} = igrp;
    end
end

% check sig orientation
[nsamps, nch] = size(sig);
if nsamps < nch
    sig = sig';
    [nsamps, nch] = size(sig);
end

% truncate sig to a multiple of fs * winstep
sig = sig(1 : (length(sig) - mod(length(sig), fs * winstep)), :);

% pad the sig signal so that the first bin starts at time 0
sig = [sig(1 : round(fs * (window - winstep) / 2), :); sig;...
    sig((end + 1 - round(fs * (window - winstep) / 2)) : end, :)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for igrp = 1 : ngrp
    [s(:, :, igrp), tstamps, freq] = mtspecgramc(sig(:, ch{igrp}),...
        [window, winstep], mtspec_params);
end

% adjust time axis
tstamps = tstamps - (window - winstep) / 2;

% adjust the fequency domain to the target frequencies by finding the
% closes indices to freq
if ~isempty(ftarget)
    fidx = zeros(1, length(ftarget)); % find closest indices in f
    for ifreq = 1 : length(ftarget)
        [fdev(ifreq), fidx(ifreq)] = min(abs(freq - ftarget(ifreq)));
    end
    spec.info.freqOrig = freq;
    freq = ftarget;
    s = s(:, fidx, :);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate power in specific bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% should use decibels
% maybe standardize with calibration data from accusleep
% increase freq resolution in lower bands by recalculating spec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% organize struct
spec.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
spec.info.tapers = mtspec_params.tapers;
spec.info.window = window;
spec.info.winstep = winstep;
spec.info.pad = mtspec_params.pad;
spec.s = s;
spec.freq = freq;
spec.tstamps = tstamps;

if saveVar
    save(specfile, 'spec')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics    
    plot_spec(spec, logfreq, basepath)  
end

end

% EOF

