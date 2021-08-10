function [psdStates, faxis, emgRMS] = psd_states(varargin)

% calculates the lfp /eeg psd for each state epoch, averages and
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
%   emg             numeric. emg data (1 x n)
%   eeg             numeric. eeg data (1 x n)
%   labels          numeric. 
%   fs              numeric. sampling frequency {1250}. 
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
% 26 jul 21 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'eeg', [], @isnumeric);
addOptional(p, 'emg', [], @isnumeric);
addOptional(p, 'labels', [], @isnumeric);
addOptional(p, 'fs', 1250, @isnumeric);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
eeg             = p.Results.eeg;
emg             = p.Results.emg;
labels          = p.Results.labels;
fs              = p.Results.fs;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% manually alternate between accuSleep signal (60 Hz low-pass) and raw data
dataMode = 'as';
if strcmp(dataMode, 'as')
    newFs = 128;
    maxF = 50;
elseif strcmp(dataMode, 'raw')
    newFs = 1250;
    maxF = 250;
end

% subsample signals
if fs ~= newFs
    eeg = standardizeSR(eeg, fs, newFs);
    emg = standardizeSR(emg, fs, newFs);
end

% state params
[cfg_colors, cfg_names, ~] = as_loadConfig([]);
nstates = length(cfg_names);
sstates = 1 : nstates - 1;      % selected states (ignore bin)
epochLen = 1;                   % period of labels [s]    
minEpDur = 10;                   % minimum segement duration

% fft params
win = hann(2 ^ (nextpow2(2 * newFs) - 1));
noverlap = floor(0.25 * newFs);
% frequencies for psd estimate. note that both the slowest frequency and
% the frequency resolution is determined by 1 / epoch length. For example,
% to estimate frequencies in a resolution of 0.2 Hz, the minimum epoch
% duration must be 5 seconds. however, using a hamming window to smooth the
% psd also reduces the the frequency resolution. further, we omit the first
% and last bin of an epoch to assure no contamination from other states.
% thus the minEpDur was set to twice the theoretical minimum (10 s).
faxis = [0.2 : 0.2 : maxF];       
    
% convert labels to state epochs. 
for istate = 1 : nstates
    binaryVec = zeros(length(labels), 1);
    binaryVec(labels == istate) = 1;
    stateEpisodes = binary2epochs('vec', binaryVec, 'minDur', [], 'maxDur', [],...
        'interDur', [], 'exclude', false); % these are given as indices and are equivalent to seconds
    stateEpochs{istate} = stateEpisodes * epochLen;
end
epLen = cellfun(@(x) (diff(x')'), stateEpochs, 'UniformOutput', false);
nepochs = cellfun(@length, epLen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc power for each epoch separatly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psdStates = zeros(nstates - 1, length(faxis));
for istate = sstates
    for iepoch = 1 : length(epLen{istate})
        
        % use only epochs longer than 
        if epLen{istate}(iepoch) < minEpDur 
            emgRMS{istate}(iepoch) = NaN;
            nepochs(istate) = nepochs(istate) - 1;
            continue
        end
        
        % idx to epoch signal w/o first and last bin
        dataIdx = (stateEpochs{istate}(iepoch, 1) + 1) * newFs :...
            (stateEpochs{istate}(iepoch, 2) - 1) * newFs - 1;
        
        % calc power
        [pow, ~] = pwelch(eeg(dataIdx), win, noverlap, faxis, newFs);
        psdStates(istate, :) = psdStates(istate, :) + pow;                       
               
        % emg rms 
        if ~isempty(emg)
            emgRMS{istate}(iepoch) = rms(emg(dataIdx));
        end
    end
    psdStates(istate, :) = psdStates(istate, :) / nepochs(istate);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot spectral power per state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use 10*log10(psdStates) for [dB]

if graphics
    xLimit = [0 maxF];
    fh = figure;
    % raw psd
    sb1 = subplot(1, 2, 1);
    ph = plot(faxis, psdStates, 'LineWidth', 3);
    set(ph, {'color'}, cfg_colors(sstates))
    xlim(xLimit)
    xlabel('Frequency [Hz]')
    ylabel('PSD [mV^2/Hz]')
    if strcmp(dataMode, 'raw')
        set(gca, 'YScale', 'log')
    end
    % norm psd
    sb2 = subplot(1, 2, 2);
    ph = plot(faxis, psdStates ./ sum(psdStates, 2), 'LineWidth', 3);
    set(ph, {'color'}, cfg_colors(sstates))
    xlim(xLimit)
    xlabel('Frequency [Hz]')
    ylabel('norm PSD')
    if strcmp(dataMode, 'raw')
        set(gca, 'YScale', 'log')
    end
end
