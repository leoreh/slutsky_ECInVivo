function [psdBins, faxis] = psd_timebins(varargin)

% calculates the psd of a signal with the pwelch method in specific time
% bins by averaging. 
%
% INPUT:
%   sig             numeric. eeg data (1 x n)
%   fs              numeric. sampling frequency {1250} of the signals.
%   faxis           numeric. frequencies of psd estimate {[0.5 : 0.2 : 100]}
%   winCalc         cell of n x 2 mats. the psd will be calculated for each
%                   cell separately according to the time windows specified
%                   in each cell [sec]
%   graphics        logical. plot figure {true}
% 
% OUTPUT
%   psdBins         raw averaged psd for each state (mat nstates x faxis)
%   faxis           frequencies of psd estimate 
%
% DEPENDENCIES
% 
% TO DO LIST
%
% 26 jul 21 LH      updates:
% 05 jan 21         removed downsampling
% 12 jan 22         winCalc. removed state dependancy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'sig', [], @isnumeric);
addOptional(p, 'fs', 1250, @isnumeric);
addOptional(p, 'faxis', [], @isnumeric);
addOptional(p, 'winCalc', []);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
sig             = p.Results.sig;
fs              = p.Results.fs;
faxis           = p.Results.faxis;
winCalc         = p.Results.winCalc;
graphics        = p.Results.graphics;

if isempty(faxis)
    faxis = [0.5 : 0.2 : 100];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assert time bins
if isempty(winCalc)
    winCalc = {[0 Inf]};
end
if ~iscell(winCalc)
    if size(winCalc, 2) ~= 2
        error('check the input winCalc')
    end
    winCalc = {winCalc};
end

maxsec = floor(length(sig) / fs);
if winCalc{end}(end) > maxsec
    winCalc{end}(end) = maxsec;
end

nwin = size(winCalc, 2);
nbins = cellfun(@(x) size(x, 1), winCalc, 'uni', true);
for iwin = 1 : nwin
    if nbins(iwin) < 1
        fprintf('\nWARNING: some timebins are empty\n')
        winCalc{iwin} = [];
    end
end

% fft params
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc power for each epoch separatly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psdBins = zeros(nwin, length(faxis));
for iwin = 1 : nwin
    for ibin = 1 : nbins(iwin)
               
        % idx to signal w/o first and last bin
        dataIdx = (winCalc{iwin}(ibin, 1) + 1) * fs :...
            (winCalc{iwin}(ibin, 2) - 1) * fs - 1;
        
        % calc power and sum across bins
        [pow, ~] = pwelch(sig(dataIdx), win, noverlap, faxis, fs);
        psdBins(iwin, :) = psdBins(iwin, :) + pow;                       
               
    end
    % average power
    psdBins(iwin, :) = psdBins(iwin, :) / nbins(iwin);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot spectral power per window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use 10*log10(psdStates) for [dB]

if graphics
    xLimit = [0 faxis(end)];
    fh = figure;
    
    % raw psd
    sb1 = subplot(1, 2, 1);
    ph = plot(faxis, psdBins, 'LineWidth', 3);
    xlim(xLimit)
    xlabel('Frequency [Hz]')
    ylabel('PSD [mV^2/Hz]')
    set(gca, 'YScale', 'log')
    
    % norm psd
    sb2 = subplot(1, 2, 2);
    ph = plot(faxis, psdBins ./ sum(psdBins, 2), 'LineWidth', 3);
    xlim(xLimit)
    xlabel('Frequency [Hz]')
    ylabel('norm PSD')
    set(gca, 'YScale', 'log')
end

end

% EOF