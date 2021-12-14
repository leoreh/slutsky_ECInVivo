function [zband, t, pband] = specBand(varargin)

% calculates power spectra in specified band. can be normalized to power in
% broadband (0-100 Hz). if band is empty then simply calculates and plots
% the broadband spectrogram.
%
% INPUT
%   sig         signal for detection
%   fs          sampling frequency {1250}.
%   band        vector. e.g. {[1 4]} for delta
%   binsize     scalar {10}. in [samples]
%   smf         smooth factor {winsize}. empty means no smoothing. zero
%               means smf = binsize
%   normband    logical. normalize to broadband (1) or not {0}
%   logaxis     logical. ploy y axis (f) on logscale {false}
%   graphics    logical. plot figure {1}
%
% OUTPUT
%   zband       power spectra in band
%   t           bin centers over which zband is calculated
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

p = inputParser;
addParameter(p, 'sig', [], @isnumeric)
addParameter(p, 'fs', 1250, @isnumeric)
addParameter(p, 'band', [1 4], @isnumeric)
addParameter(p, 'binsize', 10, @isnumeric)
addParameter(p, 'smf', 0, @isnumeric)
addParameter(p, 'normband', true, @islogical)
addParameter(p, 'logaxis', false, @islogical)
addParameter(p, 'graphics', false, @islogical)

parse(p, varargin{:})
sig = p.Results.sig;
fs = p.Results.fs;
band = p.Results.band;
binsize = p.Results.binsize;
smf = p.Results.smf;
normband = p.Results.normband;
logaxis = p.Results.logaxis;
graphics = p.Results.graphics;

% params
freq = logspace(0, 2, 100);

if smf == 0
    smf = binsize;
end

% initialize output
zband = [];
idx = zeros(1, 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% instead of 10 s window w/ 1 s overlap, 0 overlap and smooth with 10
% points. this is so that p is congruent with bsr / iis rate
win = hann(2 ^ nextpow2(binsize));
[~, f, t, pband] = spectrogram(sig, win, 0, freq, fs, 'yaxis', 'psd');

% psd in dB
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

% smooth
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
    pband = bz_NormToRange(pband, [0 1]);
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
    if logaxis
        set(gca, 'YScale', 'log')
    end
    box off
    title('Wideband spectrogram')
end

end

% EOF

