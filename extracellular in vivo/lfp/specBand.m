function [zband, t] = specBand(varargin)

% calculates power spectra in specified band. can be normalized to power in
% broadband (0-100 Hz). if band is empty then simply calculates and plots
% the broadband spectrogram.
% 
% INPUT
%   sig         signal for detection
%   fs          sampling frequency
%   band        vector. e.g. {[1 4]} for delta 
%   binsize     scalar {10}. in [samples]
%   smf         smooth factor {winsize}. empty means no smoothing. zero
%               means smf = binsize
%   normband    logical. normalize to broadband (1) or not {0}
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
% 13 jan 20 LH.     updates:
% 20 feb 20 LH      normalize to broadband

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pband = inputParser;
addParameter(pband, 'sig', [], @isnumeric)
addParameter(pband, 'fs', 1250, @isnumeric)
addParameter(pband, 'band', [1 4], @isnumeric)
addParameter(pband, 'binsize', 10, @isnumeric)
addParameter(pband, 'smf', 0, @isnumeric)
addParameter(pband, 'normband', true, @islogical)
addParameter(pband, 'graphics', false, @islogical)

parse(pband, varargin{:})
sig = pband.Results.sig;
fs = pband.Results.fs;
band = pband.Results.band;
binsize = pband.Results.binsize;
smf = pband.Results.smf;
normband = pband.Results.normband;
graphics = pband.Results.graphics;

% params
freq = logspace(0, 2, 100);

if smf == 0
    smf = binsize;
end

% initialize output
zband = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% instead of 10 s window w/ 1 s overlap, 0 overlap and smooth with 10
% points. this is so that p is congruent with bsr / iis rate
win = hann(2 ^ nextpow2(binsize));
[pband, f, t, ~] = spectrogram(sig, win, 0, freq, fs, 'yaxis', 'psd');

% psd
pband = 10 * log10(abs(pband));

% normalize to broadband
    broad = sum(pband, 1);
    [~, idx(1)] = min(abs(f - band(1)));
    [~, idx(2)] = min(abs(f - band(2)));
    narrow = sum(pband(idx(1) : idx(2), :), 1);
if normband
    zband = bz_NormToRange(zscore(narrow ./ broad), [0 1]);
%     zband = narrow ./ broad;
else
    zband = narrow;
end
if ~isempty(smf)
    zband = movmean(zband, smf);
    pband = movmean(pband, smf);
end
zband = zband(:);

% % smooth and z-score
% if ~isempty(band)
%     zband = sum(pband, 1);
%     zband = movmean(zband, smf);
%     if znorm
%         zband = bz_NormToRange(zscore(zband), [0 1]);
%     end
% else
%     zband = pband;
% end
% zband = zband(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ~same as Matlab default plot spectrogram
if graphics    
    surf(t / 60, f, pband, 'EdgeColor', 'none');
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

