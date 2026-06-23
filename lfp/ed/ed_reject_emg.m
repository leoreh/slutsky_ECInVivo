function ed = ed_reject_emg(ed, emg, fs, varargin)
% ED_REJECT_EMG Flag epileptiform discharges that coincide with high EMG.
%
%   ed = ED_REJECT_EMG(ed, emg, fs, varargin)
%
%   SUMMARY:
%       Scores each detected discharge for EMG (muscle) contamination and
%       records the score plus a pass/fail mask. Unlike the legacy
%       IED.reject_emg, this does NOT flip an accepted flag; it writes the
%       per-event feature and idxQA.emg so the wrapper can seed acceptance
%       and the GUI can triage. Ports both legacy methods.
%
%   INPUTS:
%       ed          - (Struct) After ed_detect (requires .pos, .info.fs).
%       emg         - (Vec) EMG trace, same fs/length as the detection signal
%                           (used by method 'zscore'; may be [] for 'emg_rms').
%       fs          - (Num) Sampling frequency [Hz].
%       varargin    - Parameter/Value pairs:
%           'method'      - (Char) 'zscore' (default) | 'emg_rms'.
%           'emgRms'      - (Vec)  1-Hz log-RMS vector (e.g. sSig.emg_rms),
%                                  required for method 'emg_rms'.
%           'winSec'      - (Num)  Window around each event to measure EMG [s]. {0.05}
%           'baselineSec' - (Num)  Moving baseline for the z-score [s]. {5}
%           'thrZ'        - (Num)  Pass if EMG z-score <= this. {3}
%           'thrRms'      - (Num)  Pass if emg_rms <= this. {75th percentile}
%
%   OUTPUTS:
%       ed          - (Struct) With added fields:
%           .emgZ      - (N x 1) EMG z-score per event (method 'zscore').
%           .emgRms    - (N x 1) emg_rms sampled per event (method 'emg_rms').
%           .idxQA.emg - (N x 1) logical, true = passes the EMG criterion.
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Created: 22 Jun 2026

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'ed', @isstruct);
addRequired(p, 'emg', @(x) isempty(x) || isnumeric(x));
addRequired(p, 'fs', @isnumeric);
addParameter(p, 'method', 'zscore', @(x) any(strcmpi(x, {'zscore', 'emg_rms'})));
addParameter(p, 'emgRms', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'winSec', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'baselineSec', 5, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'thrZ', 3, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'thrRms', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));

parse(p, ed, emg, fs, varargin{:});
ed     = p.Results.ed;
emg    = p.Results.emg;
fs     = p.Results.fs;
method = lower(p.Results.method);

%% ========================================================================
%  SCORE EVENTS
%  ========================================================================

pos     = ed.pos(:);
nEvents = numel(pos);
ed.emgZ   = nan(nEvents, 1);
ed.emgRms = nan(nEvents, 1);

switch method
    case 'zscore'
        if isempty(emg)
            error('ed_reject_emg:noEmg', 'Pass an EMG trace for method ''zscore''.')
        end
        emg = emg(:);
        winSamp  = max(1, round(fs * p.Results.winSec));
        baseSamp = max(1, round(fs * p.Results.baselineSec));
        emgAbs = abs(emg);
        emgMu  = movmean(emgAbs, baseSamp);
        emgSd  = movstd(emgAbs, baseSamp);
        emgSd(emgSd == 0) = eps;

        peakEmg = zeros(nEvents, 1);
        for iEv = 1:nEvents
            seg = emgAbs(max(1, pos(iEv) - winSamp) : min(numel(emgAbs), pos(iEv) + winSamp));
            peakEmg(iEv) = max(seg);
        end
        score   = (peakEmg - emgMu(pos)) ./ emgSd(pos);
        ed.emgZ = score;
        pass    = score <= p.Results.thrZ;
        thrMsg  = sprintf('thrZ=%g', p.Results.thrZ);

    case 'emg_rms'
        emgRms = p.Results.emgRms(:);
        if isempty(emgRms)
            error('ed_reject_emg:noRms', 'Pass emgRms (e.g. sSig.emg_rms) for method ''emg_rms''.')
        end
        if isempty(p.Results.thrRms)
            thrRms = prctile(emgRms, 75);
        else
            thrRms = p.Results.thrRms;
        end
        eventSec  = (pos - 1) / fs;
        score     = interp1((1:numel(emgRms))', emgRms, eventSec, 'nearest', 'extrap');
        ed.emgRms = score;
        pass      = score <= thrRms;
        thrMsg    = sprintf('thrRms=%.3f', thrRms);
end

%% ========================================================================
%  RECORD
%  ========================================================================

if ~isfield(ed, 'idxQA') || ~isstruct(ed.idxQA)
    ed.idxQA = struct();
end
ed.idxQA.emg = pass(:);
ed.info.emgMethod = method;
ed.info.emgThr = thrMsg;

fprintf('ed_reject_emg (%s): %d / %d events flagged high-EMG (%s)\n', ...
    method, sum(~pass), nEvents, thrMsg);

end     % EOF
