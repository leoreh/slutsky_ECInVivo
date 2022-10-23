function spec = calc_spec(varargin)

% calculates the multitaper spectrogram via chronux with refinements from
% accusleep. data can be loaded from .lfp or provided as input.
% multichannel support includes averaging channels (e.g. electrodes from
% the same tetrode) or avergaing spectrograms across channels and/or
% calculating the spectrogram separately for groups of channels. these are
% determined by the input 'ch'
% see wonderful explanation of multitaper estimates here:
% https://www.youtube.com/watch?v=6qTD7qtHius
%
% INPUT
%   basepath    char. fullpath to recording folder {pwd}
%   sig         lfp signal. if matrix assumes columns are samples and rows
%               are channels. if empty will load from (basename).lfp
%               according to ch
%   ch          cell of vecs depicting groups of channels (indices to the
%               rows of sig) whose spectrogram should be averaged. this is
%               done in mtspectrumc.m. if sig is loaded from binary then
%               averaging will be done on the lfp traces. for example, this
%               can be used for tetrodes (i.e. ch = spkgrp) and/or to
%               compare eeg w/ lfp signals
%   fs          sampling frequency {1250}.
%   ftarget     numeric. requested frequency range and resolution. this 
%               determines the range passed to chronux but the actual
%               resolution depends (1) on the time resolution, i.e. winstep
%               and fs, and (2) the degree of zero padding chunks of the
%               signal ('pad'). the output of chronux is interpolated by
%               nearest neighbors to the requested frequencies
%   padfft      numeric. zero pad chunks of the data to the next pow2
%               when applying fft. 0 means yes -1 means no. see mtspecgramc.m
%   winstep     numeric. determines the time resolution of the spectrogram.
%               for accusleep should be equal to epoch length. {1} [sec]
%   logfreq     logical. plot freq axis on logscale {false}. requires
%               that ftarget be log-spaced
%   force       logical. re-analyze even if spec file exists
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
%       # find a way to set request freqs in a log scale (done)
%       # normalize spectrogram (~done)
%       # separate graphics to stand alone (done)
%       # calc power in bands (e.g. delta / theta) (done)
%       # overcome tradeoff between freq and time resolution. perhaps by
%       calculating spec separately for low frequencies or by using wavelet
%
% 13 jan 20 LH      updates:
% 20 feb 20         normalize to broadband
% 25 feb 22         adapted mtspecgramc
% 18 apr 22         multichannel support
% 26 apr 22         bands
% 29 may 22         band through sum instead of mean

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
addParameter(p, 'winstep', 5, @isnumeric)
addParameter(p, 'logfreq', false, @islogical)
addParameter(p, 'force', false, @islogical)
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
force           = p.Results.force;
saveVar         = p.Results.saveVar;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% files
[~, basename] = fileparts(basepath);
specfile = fullfile(basepath, [basename, '.spec.mat']);
sessionfile = fullfile(basepath, [basename, '.session.mat']);

% load if exists
if exist(specfile, 'file') && ~force
    load(specfile, 'spec')
    return
end

% session params
if exist(sessionfile, 'file')
    load(sessionfile, 'session')
    nchans = session.extracellular.nChannels;
end

% prep frequencies
if isempty(ftarget)
    if logfreq
        ftarget = logspace(log10(0.5), 2, 200);
    else
        ftarget = linspace(0.5, 100, 200);
    end
end
frange = [ftarget(1), ftarget(end)];
if frange(2) > fs / 2
    error('requested max frequency greater than nyquist')
end

% chronux params
window = max([5, winstep]);
mtspec_params.pad = padfft;
mtspec_params.Fs = fs;
mtspec_params.fpass = frange;
mtspec_params.tapers = [ceil(window / 2), ceil(window / 2) * 2 - 1];
mtspec_params.trialave = 1;
mtspec_params.err = [0, 0];

% organize channel groups
spec.info.ch = ch;
if isempty(ch)
    ch = {1};
end
ngrp = length(ch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load and average within channel groups
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

clear s
for igrp = 1 : ngrp
    [s(:, :, igrp), tstamps, freq] = mtspecgramc(sig(:, ch{igrp}),...
        [window, winstep], mtspec_params);
end

% adjust time axis
tstamps = tstamps - (window - winstep) / 2;

% adjust the fequency domain to the target frequencies by finding the
% closest indices to freq
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

% bands taken from Boyce et al., Science, 2016
bandNames = ["broad", "swa", "delta", "theta", "alpha", "beta", "glow", "ghigh"];
bandFreqs = [0.5, 100; 0.5, 1; 1, 4; 4, 10; 10, 14; 15, 30; 30, 60; 60, 100];

% convert to dB (chronux output is power not magnitude)
powdb = 10 * log10(s);

% calc power in band. db is a 3d array of freqBand x time x channel
for iband = 1 : length(bandFreqs)
    bandIdx = InIntervals(freq, bandFreqs(iband, :));
    spec.bands.db(iband, :, :) = squeeze(sum(powdb(:, bandIdx, :), 2));
end
spec.bands.bandNames = bandNames;
spec.bands.bandFreqs = bandFreqs;

% standardization can be done by deviding with broadband or perhaps with
% the calibration data from accusleep. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% organize struct
spec.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
spec.info.tapers = mtspec_params.tapers;
spec.info.window = window;
spec.info.winstep = winstep;
spec.info.pad = mtspec_params.pad;
spec.info.fdev = fdev;
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
    % manual selections
    ch = 1;
    dataType = 'raw';   % can be raw / norm / none
    
    clear axh
    fh = figure;
    th = tiledlayout(2, 1, 'TileSpacing', 'Compact');
    axh(1) = nexttile;
    plot_spec(spec, 'ch', ch, 'logfreq', logfreq, 'saveFig', true,...
        'axh', axh)

    axh(2) = nexttile;
    for iband = 1 : length(spec.bands.bandFreqs)
        lgd{iband} = sprintf('%s [%d-%d Hz]',...
            spec.bands.bandNames{iband}, floor(spec.bands.bandFreqs(iband, 1)),...
            spec.bands.bandFreqs(iband, 2));
    end
    yval = movmean(squeeze(spec.bands.db(:, :, ch)), 100,  2);
    switch dataType
        case 'raw'
            plot(spec.tstamps / 60 / 60,...
                yval(2 : end, :), 'LineWidth', 2)
            ylabel('Spectral Power [mV2/Hz]')
        case 'norm'
            plot(spec.tstamps / 60 / 60,...
                yval(2 : end, :) ./ yval(1, :), 'LineWidth', 2)
            ylabel('Norm. Spectral Power [dB]')
        case 'none'
            return
    end
    legend(lgd{2 : end})
    xlabel('Time [hr]')
    linkaxes(axh, 'x')
end

end

% EOF

