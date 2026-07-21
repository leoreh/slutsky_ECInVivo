function otl = lfp_artifacts(sig, fs, varargin)

% Detects gross amplitude artifacts (movement steps) in a raw trace.
%
% A severe movement artifact is a step: the voltage shifts tens of times the
% normal envelope within a sample or two, then relaxes back over hundreds of ms.
% It is found here by absolute deviation from the median, in units of a robust
% (MAD-based) spread, so one constant transfers across recordings whose absolute
% scale differs - on the MCU cohort the same thrFactor lands anywhere from ~0.9
% to ~1.8 mV. That level is far above any physiological event, which is what
% lets the mask be applied blind: it cannot remove ripples or sharp waves.
%
% A derivative criterion (|diff|) was tested and rejected. Sample-to-sample jumps
% are dominated by ordinary fast LFP, so at any useful sensitivity it flags MORE
% of a clean recording than of a contaminated one. What separates an artifact is
% its size and its slow recovery; size alone is enough, and costs no parameters.
%
% TWO RESOLUTIONS, BECAUSE THE CONSUMERS DIFFER
% - .boolean / .bouts  are PADDED: they cover the deflection AND the recovery
%   that follows it. Use these to keep artifacts out of a variance / baseline
%   estimate, where the decaying tail still contributes energy long after the
%   trace has dropped back under threshold (see ripp_sigPrep).
% - .epoch / .epochBouts are UNPADDED and epoch-wise: an epoch is flagged when
%   it actually contains a deflection. Use these for spectrogram-derived work,
%   where the damage is the transient itself corrupting that window's spectrum,
%   and where padding would blank neighbouring epochs for no reason.
%
% INPUTS
% - sig             <vec> raw signal, e.g. LFP in uV.
% - fs              <num> sampling frequency [Hz].
% - thrFactor       <num>(opt) deviation from the median, in robust SDs, above
%                   which a sample is an artifact. {8}
% - pad             <num>(opt) seconds added either side of each deflection for
%                   .boolean / .bouts, covering the recovery. {0.5}
% - epochLen        <num>(opt) epoch length [s] for .epoch / .epochBouts, to
%                   match a spectrogram / label vector. {1}
% - mask            <log>(opt) restrict the robust statistics to these samples
%                   (e.g. NREM only). Flags are returned for the whole signal
%                   either way. {[] = every sample}
% - flgPlot         <log>(opt) plot the mean artifact waveform + examples, so
%                   what was detected can be confirmed by eye. {false}
%
% OUTPUT
% - otl             <struct> .boolean .bouts .epoch .epochBouts .thr .frac .info
%
% DEPENDENCIES
% - binary2bouts
%
% SEE ALSO
% - get_otl (feature-space outliers across bouts; a different job)
%
% HISTORY
% - 260720          created as get_otlAmp, for ripple baseline estimation: in one
%                   CAG session 0.8% of NREM held 94% of the detection signal's
%                   energy, inflating its scale 4x and suppressing detection
%                   ~28-fold.
% - 260720b         renamed lfp_artifacts and made the single artifact entry
%                   point (get_otlSpec retired). Split padded (baseline) from
%                   unpadded epoch (spectrogram) output; added flgPlot.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'sig', @isnumeric);
addRequired(p, 'fs', @isnumeric);
addParameter(p, 'thrFactor', 8, @isnumeric);
addParameter(p, 'pad', 0.5, @isnumeric);
addParameter(p, 'epochLen', 1, @isnumeric);
addParameter(p, 'mask', [], @(x) isempty(x) || islogical(x));
addParameter(p, 'flgPlot', false, @islogical);
parse(p, sig, fs, varargin{:});

thrFactor = p.Results.thrFactor;
pad       = p.Results.pad;
epochLen  = p.Results.epochLen;
mask      = p.Results.mask;
flgPlot   = p.Results.flgPlot;

sig = sig(:);
nSamp = numel(sig);
if isempty(mask), mask = true(nSamp, 1); else, mask = mask(:); end

otl = struct('boolean', false(nSamp, 1), 'bouts', zeros(0, 2), ...
    'epoch', false(0, 1), 'epochBouts', zeros(0, 2), 'thr', NaN, 'frac', 0);
otl.info = struct('thrFactor', thrFactor, 'pad', pad, 'epochLen', epochLen, ...
    'loc', NaN, 'scl', NaN, 'nEvents', 0, 'calcTime', datetime("now"));

if ~any(mask), return; end

%% ========================================================================
%  DETECT
%  ========================================================================
% median / MAD rather than mean / SD: the statistic must not be moved by the
% very transients it is meant to find.

loc = median(sig(mask), 'omitnan');
scl = 1.4826 * median(abs(sig(mask) - loc), 'omitnan');
if ~isfinite(scl) || scl == 0, return; end     % flat trace -> nothing to flag

thr = thrFactor * scl;
hit = sig > loc + thr | sig < loc - thr;       % no abs(): avoids a full copy
if ~any(hit)
    otl.epoch = false(floor(nSamp / round(epochLen * fs)), 1);
    otl.thr = thr; otl.info.loc = loc; otl.info.scl = scl;
    return
end

%% ========================================================================
%  TWO RESOLUTIONS
%  ========================================================================
% Epoch flags come from the raw deflections: a spectrogram window is ruined by
% the transient inside it, and padding would condemn its neighbours as well.
% The sample mask is padded, because a variance estimate is still polluted by
% the recovery tail after the trace falls back under threshold.

epochSamp = max(1, round(epochLen * fs));
nEpoch = floor(nSamp / epochSamp);
if nEpoch > 0
    hitCrop = hit(1 : nEpoch * epochSamp);
    otl.epoch = any(reshape(hitCrop, epochSamp, nEpoch), 1)';
    otl.epochBouts = boutsOf(otl.epoch);
end

hitPad = hit;
if pad > 0
    win = 2 * round(pad * fs) + 1;
    hitPad = movmax(double(hit), win) > 0;
end

otl.boolean      = hitPad;
otl.bouts        = boutsOf(hitPad) / fs;       % seconds
otl.thr          = thr;
otl.frac         = mean(hitPad(mask));
otl.info.loc     = loc;
otl.info.scl     = scl;
otl.info.nEvents = size(otl.bouts, 1);

%% ========================================================================
%  GRAPHICS
%  ========================================================================

if flgPlot
    % the PADDED events, so one waveform is one artifact: a burst of threshold
    % crossings is a single event, and n matches info.nEvents
    plotArtifacts(sig, fs, hitPad, loc, thr);
end

end     % MAIN


% =========================================================================
%  LOCALS
% =========================================================================

function b = boutsOf(vec)
% contiguous true runs of a logical vector, as [start stop] indices
b = zeros(0, 2);
if ~any(vec), return; end
b = binary2bouts('vec', vec(:), 'minDur', [], 'maxDur', [], ...
    'interDur', [], 'exclude', false, 'flgPrnt', false);
end


function plotArtifacts(sig, fs, hit, loc, thr)
% the mean artifact waveform, aligned on each event's largest deviation, with
% single events behind it. This is the figure to judge whether what was caught
% is really a movement step and not something physiological.

evt = binary2bouts('vec', hit, 'minDur', [], 'maxDur', [], 'interDur', [], ...
    'exclude', false, 'flgPrnt', false);
if isempty(evt), return; end

winSec = 1;
halfW  = round(winSec * fs);
nEvt   = size(evt, 1);

pk = nan(nEvt, 1);
for iEvt = 1 : nEvt
    seg = sig(evt(iEvt, 1) : evt(iEvt, 2));
    [~, iMax] = max(abs(seg - loc));
    pk(iEvt) = evt(iEvt, 1) + iMax - 1;
end
pk = pk(pk > halfW & pk < numel(sig) - halfW);
if isempty(pk), return; end

wv = zeros(numel(pk), 2 * halfW + 1);
for iPk = 1 : numel(pk)
    wv(iPk, :) = sig(pk(iPk) - halfW : pk(iPk) + halfW) - loc;
end
wv = wv .* sign(wv(:, halfW + 1));             % align polarity on the peak
tms = (-halfW : halfW) / fs * 1000;

fh = figure('Color', 'w');
tiledlayout(fh, 1, 2, 'TileSpacing', 'compact');

nexttile
plot(tms, wv(1 : min(end, 200), :)', 'Color', [0.7 0.7 0.7 0.3]);
hold on
plot(tms, mean(wv, 1), 'k', 'LineWidth', 2);
yline([thr, -thr], '--r');
xlabel('Time from peak (ms)'); ylabel('Deviation from median');
title(sprintf('Artifact waveform (n = %d)', numel(pk)));
axis tight

nexttile
histogram(max(abs(wv), [], 2) / thr, 40);
xlabel('Peak deviation (x threshold)'); ylabel('Count');
title('Amplitude distribution');
set(gca, 'xscale', 'log')

end     % EOF
