function [psdStates, faxis, epStats] = psd_states(varargin)

% calculates the lfp / eeg psd for each state epoch, averages and
% normalizes. can also calculate the emg rms for each epoch. subsamples
% signals to 128 Hz for faster computation time. this fits lfp data from
% as_prepSig where a low-pass of 60 Hz is applied. if higher frequencies
% are wanted then use raw data and subsample to nyquist (or do not
% subsample at all).

% another parameter observed in Yuval Nir's code is eeg ratio high / low
% calculated as: eegRatio{istate}(iepoch) = sum(pow(faxis > 25)) /
% sum(pow(faxis < 5));
%
% INPUT:
%   eeg             numeric. eeg data (1 x n)
%   labels          numeric. 
%   fs              numeric. sampling frequency {1250} of the signals.
%   faxis           numeric. frequencies of psd estimate {[0.2 : 0.2 : 120]}
%   sstates         numeric. idx of selected states. e.g. [1, 4, 5] will
%                   only plot wake, nrem and rem
%   graphics        logical. plot figure {true}
% 
% OUTPUT
%   psdStates       raw averaged psd for each state (mat nstates x faxis)
%   faxis           frequencies of psd estimate 
%   emgRMS          rms for each epoch, nstates x epochs (cell of nstates)
%
% DEPENDENCIES
% 
% TO DO LIST
%   allow input of ss struct or stateEpochs instead of labels
%
% 26 jul 21 LH      updates:
% 05 jan 21         removed downsampling  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'eeg', [], @isnumeric);
addOptional(p, 'labels', [], @isnumeric);
addOptional(p, 'fs', 1250, @isnumeric);
addOptional(p, 'faxis', [0.2 : 0.2 : 120], @isnumeric);
addOptional(p, 'sstates', [], @isnumeric);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
eeg             = p.Results.eeg;
labels          = p.Results.labels;
fs              = p.Results.fs;
faxis           = p.Results.faxis;
sstates         = p.Results.sstates;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fft params
if isempty(faxis)
    faxis = [0.2 : 0.2 : 120];       
end
win = hann(2 ^ (nextpow2(2 * fs) - 1));
noverlap = floor(0.25 * fs);
% frequencies for psd estimate. note that both the slowest frequency and
% the frequency resolution is determined by 1 / epoch length. For example,
% to estimate frequencies in a resolution of 0.2 Hz, the minimum epoch
% duration must be 5 seconds (minDur). however, using a hamming window to
% smooth the psd also reduces the frequency resolution. further, we omit
% the first and last bin of an epoch to assure no contamination from other
% states. this is why the minDur was set to twice the theoretical minimum
% (10 s).

% state params
cfg = as_loadConfig();
nstates = cfg.nstates;
if isempty(sstates)
    sstates = 1 : nstates;      % selected states (ignore bin)
end

% convert labels to state epochs. 
minDur = [10, 5, 5, 10, 5, 5];  % threshold of minimum epoch length
interDur = 4;                   % combine epochs separated by <= interDur

if length(minDur) == 1
    minDur = repamt(minDur, nstates, 1);
elseif length(minDur) ~= nstates
    error('minDur length is different than the number of states')
end
if length(interDur) == 1
    interDur = repmat(interDur, nstates, 1);
elseif length(interDur) ~= nstates
    error('interDur length is different than the number of states')
end

% create state epochs
for istate = 1 : nstates
    binaryVec = zeros(length(labels), 1);
    binaryVec(labels == istate) = 1;
    stateEpochs{istate} = binary2epochs('vec', binaryVec, 'minDur', minDur(istate), 'maxDur', [],...
        'interDur', interDur(istate), 'exclude', false);
end

% epoch stats
epStats.epLen = cellfun(@(x) (diff(x')'), stateEpochs, 'UniformOutput', false);
epStats.nepochs = cellfun(@length, epStats.epLen);
epStats.totDur = cellfun(@sum, epStats.epLen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc power for each epoch separatly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psdStates = zeros(length(sstates), length(faxis));
for istate = 1 : length(sstates)
    sidx = sstates(istate);
    for iepoch = 1 : epStats.nepochs(sidx)
               
        % idx to epoch signal w/o first and last bin
        dataIdx = (stateEpochs{sidx}(iepoch, 1) + 1) * fs :...
            (stateEpochs{sidx}(iepoch, 2) - 1) * fs - 1;
        
        % calc power and sum across epochs
        [pow, ~] = pwelch(eeg(dataIdx), win, noverlap, faxis, fs);
        psdStates(istate, :) = psdStates(istate, :) + pow;                       
               
    end
    % average power
    psdStates(istate, :) = psdStates(istate, :) / epStats.nepochs(sidx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot spectral power per state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use 10*log10(psdStates) for [dB]

if graphics
    xLimit = [0 faxis(end)];
    fh = figure;
    % raw psd
    sb1 = subplot(1, 2, 1);
    ph = plot(faxis, psdStates, 'LineWidth', 3);
    set(ph, {'color'}, cfg.colors(sstates))
    xlim(xLimit)
    xlabel('Frequency [Hz]')
    ylabel('PSD [mV^2/Hz]')
    set(gca, 'YScale', 'log')
    
    % norm psd
    sb2 = subplot(1, 2, 2);
    ph = plot(faxis, psdStates ./ sum(psdStates, 2), 'LineWidth', 3);
    set(ph, {'color'}, cfg.colors(sstates))
    xlim(xLimit)
    xlabel('Frequency [Hz]')
    ylabel('norm PSD')
    set(gca, 'YScale', 'log')
end
