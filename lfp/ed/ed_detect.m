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
%       Detection runs on 60-150 Hz because a discharge is defined by being
%       SHARP, not by being large. The obvious choice - a band around the
%       deflection itself, say 10-100 Hz - fails, and fails in a way worth
%       recording: in the slow band a discharge and an ordinary hippocampal
%       sharp wave are not separable at all (AUC 0.54 against curated labels),
%       so a detector built there proposes thousands of normal deflections and
%       buries the few real discharges. Above 60 Hz the same events separate
%       almost perfectly (AUC 0.99). The band also sits above 50 Hz mains and
%       is kept below the very high frequencies where EMG dominates.
%
%       Detection is polarity-blind and deliberately permissive: it thresholds
%       |filt|, so a candidate costs one row while a missed discharge is gone
%       for good. Shape is the gate's job (ed_params .posZ), not detection's.
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
%                      .peakTime .pos .bouts .fastZ .posZ .amp .dur .emg
%                      .state .accepted, partial .info).
%       aux - <struct> .edSig .fs .edCh - the prepared signal, so the caller can
%                      build the per-event maps without reloading.
%
%   DEPENDENCIES:
%       basepaths2vars, ed_methods, ed_sigLoad, ed_params, filterLFP,
%       binary2bouts, evt_boutTimes, evt_emgScore, evt_states.
%
%   HISTORY:
%       Created: 260622 (as a signal-in / events-out detector).
%       Updated: 260720 (stage 1 of the staged pipeline).
%       Updated: 260721 (detection band moved from 10-100 Hz to 60-150 Hz after
%                curated discharges showed the two are not separable below
%                ~30 Hz; the block-wise baseline and the ring features went
%                with it. See dev/ed_pipeline_rebuild.md.)

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
ed.emg = evt_emgScore(emg, [ed.peakTime - 0.02, ed.peakTime + 0.02], fs, ...
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
ed.info.thr      = met.thr;
ed.info.limDur   = met.limDur;

aux = struct('edSig', edSig, 'fs', fs, 'edCh', edCh);

end     % EOF


% =========================================================================
%  LOCAL
% =========================================================================
function edSig = sigPrep(raw, fs, met)
% Band-limit the trace and set the two scales the features are measured in.
%
% Both scales are one median absolute deviation over the whole recording, taken
% on a stride - a median and a MAD are stable under decimation, and a fixed
% stride keeps the value reproducible run to run. A recording-wide scale (not a
% moving one) is deliberate: it does not adapt away the very quiet background a
% discharge stands out from, and a robust estimator means the discharges
% themselves cannot inflate the scale that measures them.
%
% .filt is left unsmoothed. A window wide enough to steady a threshold crossing
% also halves the peak of a sharp transient, which cost a third of the confirmed
% discharges when it was tried.
raw = double(raw(:));
stride = max(1, floor(numel(raw) / 2e6));

edSig.lfp  = raw;
edSig.filt = filterLFP(raw, 'fs', fs, 'type', 'butter', 'dataOnly', true, ...
    'order', 3, 'passband', met.passband, 'graphics', false);

edSig.sclFast = robustScale(edSig.filt(1 : stride : end));
edSig.sclRaw  = robustScale(raw(1 : stride : end));
edSig.z = abs(edSig.filt) / edSig.sclFast;

end     % sigPrep


function s = robustScale(x)
% Median absolute deviation, scaled to a standard deviation; 1 if degenerate.
s = 1.4826 * median(abs(x - median(x)));
if ~isfinite(s) || s <= 0, s = 1; end

end     % robustScale
