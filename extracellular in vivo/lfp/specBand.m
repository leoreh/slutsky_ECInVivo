function [zband, tband] = specBand(varargin)

% calculates normalized z score of power spectra in specified band. if band
% is empty then simply calculates and plots the broadband spectrogram. can
% also return broadband if band not specified.
% 
% INPUT
%   sig         signal for detection
%   fs          sampling frequency
%   band        vector. e.g. {[1 4]} for delta 
%   binsize     scalar {10}. in [samples]
%   smf         smooth factor {winsize}. empty means no smoothing
%   znorm       logical. normalize z-score (1) or not {0}
%   graphics    logical. plot figure {1}
%
% OUTPUT
%   zband       power spectra in band
%
% TO DO LIST
%       # make compatible with band as matrix (multiple bands)
%
% CALLS
%       spectrogram     
%
% 13 jan 20 LH.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'sig', [], @isnumeric)
addParameter(p, 'fs', 1250, @isnumeric)
addParameter(p, 'band', [1 4], @isnumeric)
addParameter(p, 'binsize', 10, @isnumeric)
addParameter(p, 'smf', 0, @isnumeric)
addParameter(p, 'znorm', true, @islogical)
addParameter(p, 'graphics', false, @islogical)

parse(p, varargin{:})
sig = p.Results.sig;
fs = p.Results.fs;
band = p.Results.band;
binsize = p.Results.binsize;
smf = p.Results.smf;
znorm = p.Results.znorm;
graphics = p.Results.graphics;

% params
if isempty(band)
    freq = logspace(0, 2, 100);
else
    freq = linspace(band(1), band(2), 40);
end
if smf == 0
    smf = binsize;
end
zband = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

win = hann(2 ^ nextpow2(binsize));

% instead of 10 s window w/ 1 s overlap, 0 overlap and smooth with 10
% points. this is so that bs and p are congruent
[~, f, tband, pband] = spectrogram(sig, win, 0, freq, fs, 'yaxis', 'psd');

% psd
pband = 10 * log10(abs(pband));

% smooth and z-score
if ~isempty(band)
    zband = sum(pband, 1);
    zband = movmean(zband, smf);
    if znorm
        zband = bz_NormToRange(zscore(zband), [0 1]);
    end
else
    zband = pband;
end
zband = zband(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ~same as Matlab default plot spectrogram
if graphics    
    surf(tband / 60, f, pband, 'EdgeColor', 'none');
    axis xy;
    axis tight;
    view(0,90);
    origSize = get(gca, 'Position');
    colormap(jet);
    colorbar;
    ylabel('Frequency [Hz]');
    set(gca, 'Position', origSize);
    set(gca, 'TickLength', [0 0])
    if isempty(band)
        set(gca, 'YScale', 'log')
    end
    box off
    title('Wideband spectrogram')
end
end

% EOF

