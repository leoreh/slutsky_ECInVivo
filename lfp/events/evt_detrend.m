function wv = evt_detrend(wv, tstamps, met)
% EVT_DETREND Remove a per-event linear baseline from event waveforms.
%
%   wv = EVT_DETREND(wv, tstamps, met)
%
%   SUMMARY:
%       One implementation of "take the slow drift out of each snippet",
%       shared by the clustering (evt_clust) and by anything plotting an
%       average waveform.
%
%       Why an average needs it. Every event sits on its own slow deflection,
%       and those offsets are unrelated to the event. Averaging raw snippets
%       averages the offsets too, so the mean is smeared vertically and its
%       flanks sag toward whatever the drifts happened to do. Detrending first
%       leaves a mean in real signal units whose baseline is flat by
%       construction.
%
%       This does NOT normalise. Amplitude is the thing a plotted average is
%       usually reporting, so it is left alone; evt_clust normalises separately
%       and only for the shape components.
%
%       'edge' fits the line on the FLANKS - the samples beyond half the
%       window - so the event cannot tilt the baseline it is measured against.
%       Fitting over the whole window instead lets a large asymmetric event
%       drag the line, and drag it by an amount that depends on its own
%       polarity and asymmetry, which is a distortion correlated with the shape
%       being measured. This is the snipFromBinary convention.
%
%   INPUTS:
%       wv      - <mat>  [nEv x nSamp] per-event waveforms, one row per event.
%       tstamps - <vec>  [1 x nSamp] window time base (s), from evt_maps.
%       met     - <char> (opt) 'edge' | 'full' | 'none'. {'edge'}
%
%   OUTPUT:
%       wv      - <mat>  [nEv x nSamp] with each row's linear trend removed.
%                        An all-NaN row (an event at a recording edge) stays
%                        NaN rather than poisoning the fit for the rest.
%
%   HISTORY:
%       260721 promoted out of evt_clust so a plotted average and the clustering
%              cannot drift apart. See lfp/ed/dev/ed_pipeline_rebuild.md.

if nargin < 3 || isempty(met), met = 'edge'; end
met = lower(met);
if strcmp(met, 'none'), return; end

t = tstamps(:);
if strcmp(met, 'edge')
    iFit = abs(t) >= 0.5 * max(abs(t));
    if nnz(iFit) < 3, iFit = true(size(t)); end
else
    iFit = true(size(t));
end

t = t - mean(t(iFit));
A = [t(iFit), ones(nnz(iFit), 1)];
B = [t, ones(numel(t), 1)];

% a row with any non-finite sample cannot be fitted; leave it as it is
iOk = all(isfinite(wv), 2);
wv(iOk, :) = wv(iOk, :) - (B * (A \ wv(iOk, iFit)'))';

end     % EOF
