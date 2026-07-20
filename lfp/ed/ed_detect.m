function [ed, aux] = ed_detect(basepath, varargin)
% ED_DETECT Detect + characterise epileptiform discharges; writes nothing.
%
%   [ed, aux] = ED_DETECT(basepath, varargin)
%
%   SUMMARY:
%       Stage 1 of the ED pipeline (detect -> curate), mirroring ripp_detect.
%       Loads the signal, band-passes it, thresholds the result into candidates,
%       and measures every per-event feature (ed_params, plus EMG and the
%       vigilance state). It makes NO acceptance decision - .accepted is seeded
%       all-true and the gate is a separate stage over the saved struct, so one
%       detection feeds any per-mouse curation. Times are window-relative; the
%       caller shifts to absolute. Nothing is written to disk.
%
%       The detection statistic is the band-passed trace over ONE robust scale
%       for the whole recording (median absolute deviation, on a stride). Band-
%       passing is what makes a single threshold mean the same thing in every
%       vigilance state: on the raw trace the NREM delta amplitude dominates the
%       variance, so a raw threshold is really a slow-wave detector. A robust
%       scale is what keeps a discharge from raising its own threshold, which
%       the previous moving mean / SD baseline did.
%
%   INPUTS:
%       basepath - <char> session directory.
%       varargin - Parameter/Value:
%           'basename' - <char>   file stem. {folder name}
%           'met'      - <struct> one ed_methods() config. {ed_methods}
%           'win'      - <vec>    window [start end] (s). {[0 Inf]}
%           'edCh'     - <num>    explicit channel for met.chMode 'ripp'. {[]}
%           'verbose'  - <log>    print progress. {false}
%
%   OUTPUTS:
%       ed  - <struct> window-relative candidates + per-event fields (.times
%                      .peakTime .pos .bouts .amp .ampG .ampZ .hfRatio .dur
%                      .emg .state .accepted, partial .info).
%       aux - <struct> .edSig .fs .edCh - the prepared signal, so the caller can
%                      build the per-event maps without reloading.
%
%   DEPENDENCIES:
%       basepaths2vars, ed_methods, ed_sigLoad, ed_params, filterLFP,
%       binary2bouts, evt_boutTimes, evt_emgScore, evt_states.
%
%   HISTORY:
%       Created: 260622 (as a signal-in / events-out detector).
%       Updated: 260720 (stage 1 of the staged pipeline. The moving mean / SD
%                z-score, the DC-dependent amplitude gate and the hand-rolled
%                crossing merge are gone. A block-wise baseline was tried and
%                dropped: it cost 45 lines and changed no result.)

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'basepath', @ischar);
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'met', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'win', [0 Inf], @isnumeric);
addParameter(p, 'edCh', [], @isnumeric);
addParameter(p, 'verbose', false, @islogical);
parse(p, basepath, varargin{:});

met     = p.Results.met;
win     = p.Results.win;
edCh    = p.Results.edCh;
verbose = p.Results.verbose;

basename = p.Results.basename;
if isempty(basename), [~, basename] = fileparts(basepath); end
if isempty(met), met = ed_methods('default'); end

%% ========================================================================
%  SIGNAL + CONTEXT
%  ========================================================================
% session may be absent on a minimal (sleep_sig only) layout, so it is guarded;
% the 'eeg' channel mode needs none of it

v = basepaths2vars('basepaths', {basepath}, ...
    'vars', {'session', 'sleep_states'});
session = [];
if isfield(v, 'session'), session = v.session; end

if isinf(win(2)), winDur = Inf; else, winDur = win(2) - win(1); end
boutTimes = evt_boutTimes(v, win, winDur);

if verbose, fprintf('[ED_DETECT] %s : loading signal...\n', basename); end
[raw, emg, fs, edCh] = ed_sigLoad(basepath, 'basename', basename, ...
    'chMode', met.chMode, 'edCh', edCh, 'win', win, 'session', session);

edSig = sigPrep(raw, fs, met);

%% ========================================================================
%  CANDIDATES
%  ========================================================================
% One binary2bouts call replaces the old chain of pairwise crossing merges,
% uniquetol de-duplication and the twin-peak rule, each of which could delete an
% event because of an unrelated neighbour.

limSamp = round(met.limDur / 1000 * fs);
bouts = binary2bouts('vec', edSig.z > met.thr, 'minDur', limSamp(1), ...
    'maxDur', limSamp(2), 'interDur', limSamp(3));
if isempty(bouts), bouts = zeros(0, 2); end

% peak = the largest deflection of the band-passed trace inside the candidate,
% by MAGNITUDE, so a negative-going discharge is localised as well as a positive
% one (discharges in this preparation are predominantly negative)
nEv = size(bouts, 1);
pos = zeros(nEv, 1);
for iEv = 1 : nEv
    idx = bouts(iEv, 1) : bouts(iEv, 2);
    [~, iRel] = max(abs(edSig.filt(idx)));
    pos(iEv) = bouts(iEv, 1) + iRel - 1;
end

ed = struct('pos', pos, 'bouts', bouts);
ed.info.fs = fs;

%% ========================================================================
%  PER-EVENT FEATURES  (curation is a separate stage)
%  ========================================================================
if verbose, fprintf('[ED_DETECT] %s : %d candidates\n', basename, nEv); end

ed = ed_params(edSig, ed);

% EMG over a fixed window about the peak: a discharge is a point event, so a
% 6 ms and a 40 ms one are scored against the same amount of muscle signal
ed.emg = evt_emgScore(emg, [ed.peakTime - 0.025, ed.peakTime + 0.025], fs, ...
    'baselineTimes', []);

% per-event state label; the per-bout rate table is a post-curation product
ed.state = evt_states(ed.times, ed.peakTime, boutTimes, ...
    'flgSave', false, 'flgPlot', false, 'name', 'ed', 'lbl', 'ED');

ed.accepted = true(nEv, 1);         % the gate runs downstream

%% ========================================================================
%  PROVENANCE (partial; the wrapper finalises to absolute time)
%  ========================================================================
ed.info.sigDur   = numel(raw) / fs;
ed.info.met      = met.name;
ed.info.chMode   = met.chMode;
ed.info.edCh     = edCh;
ed.info.passband = met.passband;
ed.info.hfBand   = met.hfBand;
ed.info.thr      = met.thr;
ed.info.limDur   = met.limDur;

aux = struct('edSig', edSig, 'fs', fs, 'edCh', edCh);

end     % EOF


% =========================================================================
%  LOCAL
% =========================================================================
function edSig = sigPrep(raw, fs, met)
% Band-limit the trace and normalise it by one robust scale.
%
% .filt is the discharge band and .hf a supra-physiological band no discharge
% reaches, which ed_params turns into the .hfRatio contamination metric. Both
% are kept unsmoothed: a window wide enough to steady a threshold crossing also
% halves the peak of a sharp discharge, which cost a third of the confirmed
% discharges when it was tried.
raw = double(raw(:));

% keep the HF band below Nyquist; a lower-rate session simply gets a narrower
% one rather than a filter-design error
hfBand = met.hfBand;
hfBand(2) = min(hfBand(2), 0.95 * fs / 2);

edSig.lfp  = raw;
edSig.filt = filterLFP(raw, 'fs', fs, 'type', 'butter', 'dataOnly', true, ...
    'order', 3, 'passband', met.passband, 'graphics', false);
if hfBand(1) < hfBand(2)
    edSig.hf = filterLFP(raw, 'fs', fs, 'type', 'butter', 'dataOnly', true, ...
        'order', 3, 'passband', hfBand, 'graphics', false);
else
    edSig.hf = nan(size(raw));
end

% one scale for the recording, on a stride - a median and a MAD are stable under
% decimation, and a fixed stride keeps the value reproducible run to run
sub = edSig.filt(1 : max(1, floor(numel(raw) / 2e6)) : end);
scl = 1.4826 * median(abs(sub - median(sub)));
if ~isfinite(scl) || scl <= 0, scl = 1; end

edSig.z = abs(edSig.filt) / scl;

end     % sigPrep
