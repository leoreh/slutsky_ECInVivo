function [s, tstamps, freq] = calc_spec(varargin)

% creates a spectrogram with the multitaper spectrogram by chronux, with
% refinements from accusleep. should also calculate the power in specific
% bands (delta, theta) during specified time windews. should normalize
%
% INPUT
%   basepath    char. fullpath to recording folder {pwd}
%   sig         signal for detection
%   fs          sampling frequency {1250}.
%   ftarget     numeric. target frequency range and resolution. this can 
%               be used to control the frequency axis of the spectrogram
%               which depends on (1) the degree of zero padding chunks of
%               the signal for the fft ('pad') and (2) the time
%               resolution (frequency binsize = 1 / window). {[]}. if empty
%               than the freuqency range will be [0 120] and the resolution
%               will be determined by window and pad.
%   padftt      numeric. zero pad chunks of the data to the next pow2
%               when applying fft. 0 means yes -1 means no. see mtspecgramc.m
%   winstep     numeric. determines the time resolution of the spectrogram.
%               for accusleep should be equal to epoch length. {1} [sec]
%   logfreq     logical. ploy y axis (freq) on logscale {false}
%   graphics    logical. plot figure {false}
%   saveVar     logical. organize and save struct {true}
%
% CALLS
%   mtspecgramc
% 
% TO DO LIST
%       # find a way to set the frequency resolution in a log scale
%       # normalize spectrogram
%       # separate graphics to stand alone
%
% 13 jan 20 LH      updates:
% 20 feb 20         normalize to broadband
% 25 feb 22         adapted mtspecgramc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar)
addParameter(p, 'sig', [], @isnumeric)
addParameter(p, 'fs', 1250, @isnumeric)
addParameter(p, 'ftarget', [], @isnumeric)
addParameter(p, 'padftt', 0, @isnumeric)
addParameter(p, 'winstep', 1, @isnumeric)
addParameter(p, 'logfreq', false, @islogical)
addParameter(p, 'graphics', false, @islogical)
addParameter(p, 'saveVar', true, @islogical)

parse(p, varargin{:})
basepath        = p.Results.basepath;
sig             = p.Results.sig;
fs              = p.Results.fs;
ftarget         = p.Results.ftarget;
padftt          = p.Results.padftt;
winstep         = p.Results.winstep;
logfreq         = p.Results.logfreq;
graphics        = p.Results.graphics;
saveVar         = p.Results.saveVar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% recording
[~, basename] = fileparts(basepath);
specfile = fullfile(basepath, [basename, '.spec.mat']);

% prep frequencies
if isempty(ftarget)
    frange = [1 120];
else
    frange = [ftarget(1), ftarget(end)];
    if frange(2) > fs / 2
        error('requested max frequency greater than nyquist')
    end
end

% mtspecgramc params
window = max([5, winstep]);
mtspec_params.pad = padftt;
mtspec_params.Fs = fs;
mtspec_params.fpass = frange;
mtspec_params.tapers = [3 5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep sig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% should add option to load data from lfp binary and average multiple
% channels. see trialave in mtspecgramc.m

% assert row
if ~isrow(sig)
    sig = sig';
end

% truncate sig to a multiple of fs * winstep
sig = sig(1 : (length(sig) - mod(length(sig), fs * winstep)));

% pad the sig signal so that the first bin starts at time 0
sig = [sig(1 : round(fs * (window - winstep) / 2)), sig,...
    sig((end + 1 - round(fs * (window - winstep) / 2)) : end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[s, tstamps, freq] = mtspecgramc(sig, [window, winstep], mtspec_params);

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
    s = s(:, fidx);
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

if saveVar
    spec.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
    spec.info.tapers = mtspec_params.tapers;
    spec.info.window = window;
    spec.info.winstep = winstep;
    spec.info.pad = mtspec_params.pad;
    spec.s = s;
    spec.freq = freq;
    spec.tstamps = tstamps;
    
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

