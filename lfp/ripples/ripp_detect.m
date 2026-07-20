function [ripp, aux] = ripp_detect(basepath, varargin)
% RIPP_DETECT Detect + characterise ripples for one method; writes nothing.
%
%   [ripp, aux] = RIPP_DETECT(basepath, varargin)
%
%   SUMMARY:
%       The detect stage of the ripple pipeline (stage 1 of detect -> curate ->
%       analyze), used by ripp_wrapper and ripp_screen. For one met it loads and
%       prepares the signal (or reuses an injected one), detects candidate events,
%       and measures every per-event feature: ripple params, the QA metrics
%       (EMG z, MUA gain), and the vigilance-state label. It does NOT
%       decide acceptance - .accepted is seeded all-true and the QA gate is a
%       separate stage (ripp_gate / ripp_curate) applied to the saved struct, so
%       the same detection feeds any per-mouse curation. Times are window-relative;
%       the caller shifts to absolute. Nothing is written to disk.
%
%   INPUTS:
%       basepath - <char> session directory.
%       varargin - Parameter/Value:
%           'met'     - <struct> one ripp_methods() config. {ripp_methods}
%           'win'     - <vec>    window [start end] (s). {[0 Inf]}
%           'v'       - <struct> pre-loaded vars (session, spikes, spktimes,
%                                sleep_states, units). {loaded}
%           'rippCh'  - <num>    explicit channel; else resolved. {[]}
%           'sig'     - <struct> injected signal bundle (.rippSig .emg .fs
%                                .rippCh) to skip load+prep; reuse across methods
%                                that share a signal config. {built here}
%           'verbose' - <log>    print progress. {false}
%
%   OUTPUTS:
%       ripp - <struct> window-relative events + per-event fields (.times
%                       .peakTime .state .accepted .emg .spkGain .ctrlTimes,
%                       ripple params, partial .info). ALL detected events;
%                       .accepted is seeded all-true (gate in ripp_curate).
%       aux  - <struct> .sig (signal bundle, reuse across methods) .spkTimes
%                       .muTimes .uType .fs .rippCh (reused by ripp_analyze).
%
%   DEPENDENCIES:
%       basepaths2vars, ripp_methods, evt_spkPrep, evt_boutTimes, ripp_pickCh,
%       ripp_sigLoad, ripp_sigPrep, ripp_noiseFloor, ripp_times, ripp_params,
%       evt_emgScore, evt_spkGain, evt_ctrlTimes, evt_states.
%
%   HISTORY:
%       260719 factored out of ripp_wrapper as the shared detection core.
%       260719b split: detect computes features only; the QA gate moved to
%               ripp_gate/ripp_curate (accepted seeded all-true here).

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'basepath', @ischar);
addParameter(p, 'met', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'win', [0 Inf], @isnumeric);
addParameter(p, 'v', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'rippCh', [], @isnumeric);
addParameter(p, 'sig', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'verbose', false, @islogical);
parse(p, basepath, varargin{:});
met     = p.Results.met;
win     = p.Results.win;
v       = p.Results.v;
rippCh  = p.Results.rippCh;
sig     = p.Results.sig;
verbose = p.Results.verbose;
[~, basename] = fileparts(basepath);

if isempty(met), met = ripp_methods('default'); end

%% ========================================================================
%  CONTEXT (session, spikes, bouts)
%  ========================================================================

if isempty(v)
    v = basepaths2vars('basepaths', {basepath}, ...
        'vars', {'session', 'spikes', 'spktimes', 'sleep_states', 'units'});
end
fsSpk = v.session.extracellular.sr;
if isinf(win(2)), sigDur = Inf; else, sigDur = win(2) - win(1); end

[spkTimes, muTimes, uType]       = evt_spkPrep(v, win, sigDur, fsSpk);
[boutTimes, vldTimes, nremTimes] = evt_boutTimes(v, win, sigDur);

%% ========================================================================
%  SIGNAL (built once; reused across methods via 'sig')
%  ========================================================================

if isempty(sig)
    if isempty(rippCh)
        rippCh = ripp_pickCh(basepath, 'basename', basename, ...
            'session', v.session, 'win', win, 'nremTimes', nremTimes, ...
            'flgForce', strcmp(met.chMode, 'best'));
    end
    if verbose
        fprintf('[RIPP_DETECT] %s : ch %s\n', basename, mat2str(rippCh));
    end
    [lfp, emg, fs] = ripp_sigLoad(basepath, 'win', win, 'session', v.session, ...
        'basename', basename, 'rippCh', rippCh, 'bit2uv', []);
    rippSig = ripp_sigPrep(lfp, fs, 'detectMet', met.detectMet, ...
        'passband', met.passband, 'zMet', met.zMet, 'nremTimes', nremTimes);
    sig = struct('rippSig', rippSig, 'emg', emg, 'fs', fs, 'rippCh', rippCh);
end
rippSig = sig.rippSig;
emg     = sig.emg;
fs      = sig.fs;
rippCh  = sig.rippCh;

%% ========================================================================
%  DETECT (threshold fixed, or calibrated to the 1/f noise floor)
%  ========================================================================

thr = met.thr;
chi = NaN;
if met.calibThr
    [thrPk, nf] = ripp_noiseFloor(rippSig.lfp, fs, met, ...
        'nremTimes', nremTimes, 'targetFP', met.targetFP);
    thr(2) = thrPk;
    thr(1) = max(0.5, thrPk - (met.thr(2) - met.thr(1)));
    chi = nf.chi;
end
ripp = ripp_times(rippSig, fs, 'thr', thr, 'limDur', met.limDur);
ripp = ripp_params(rippSig, ripp);

%% ========================================================================
%  QA METRICS + STATE LABEL (curation is a separate stage)
%  ========================================================================
% Detection computes the per-event QA metrics and the state label; it does NOT
% decide acceptance. The gate (state + metric ranges -> accepted) is applied to
% the saved struct by ripp_gate / ripp_curate, so one detection serves any
% per-mouse curation.

ripp.emg       = evt_emgScore(emg, ripp.times, fs, 'baselineTimes', nremTimes);
ripp.spkGain   = evt_spkGain(muTimes, ripp.times);
ripp.ctrlTimes = evt_ctrlTimes(ripp.times, 'vldTimes', vldTimes, 'flgPlot', false);

% per-event vigilance-state label (all events); the bout rate/density table is a
% post-curation product built by ripp_curate from the accepted mask
ripp.state = evt_states(ripp.times, ripp.peakTime, boutTimes, ...
    'flgSave', false, 'flgPlot', false, 'name', 'ripp', 'lbl', 'Ripple');

% seed accepted all-true; the gate runs downstream
ripp.accepted = true(numel(ripp.peakTime), 1);

if verbose
    fprintf('[RIPP_DETECT] %s : %d events, %d accepted (%.0f%%)\n', ...
        basename, numel(ripp.accepted), sum(ripp.accepted), ...
        100 * mean(ripp.accepted));
end

%% ========================================================================
%  PROVENANCE (partial; the wrapper finalises to absolute time)
%  ========================================================================

ripp.info.rippCh    = rippCh;
ripp.info.passband  = met.passband;
ripp.info.detectMet = met.detectMet;
ripp.info.zMet      = met.zMet;
ripp.info.thr       = thr;
ripp.info.chi       = chi;
ripp.info.met       = met.name;

aux.sig        = sig;
aux.spkTimes   = spkTimes;
aux.muTimes    = muTimes;
aux.uType      = uType;
aux.fs         = fs;
aux.rippCh     = rippCh;

end     % EOF
