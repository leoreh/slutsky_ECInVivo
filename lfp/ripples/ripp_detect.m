function [ripp, aux] = ripp_detect(basepath, varargin)
% RIPP_DETECT Detect + characterise ripples for one method; writes nothing.
%
%   [ripp, aux] = RIPP_DETECT(basepath, varargin)
%
%   SUMMARY:
%       The shared detection core of the ripple pipeline, used by both
%       ripp_wrapper (which adds spikes, saving, plotting, curation) and
%       ripp_screen (which sweeps methods). For one met it loads and prepares the
%       signal (or reuses an injected one), detects candidate events, measures
%       per-event params, LFP maps, and the QA metrics (EMG, MUA gain), labels
%       vigilance state, and sets the per-event acceptance mask. QA MARKS, it does
%       not remove: every detected event stays in the struct with an .accepted
%       flag, so a gate threshold can be re-screened from the saved metrics and
%       wake events remain inspectable. Times are window-relative; the caller
%       shifts to absolute. Nothing is written to disk.
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
%           'mapDur'  - <vec>    LFP map window [pre post] (s). {[-0.1 0.1]}
%           'verbose' - <log>    print progress. {false}
%
%   OUTPUTS:
%       ripp - <struct> window-relative events + per-event fields (.times
%                       .peakTime .state .accepted .emg .spkGain .ctrlTimes,
%                       ripple params, partial .info). ALL detected events; use
%                       .accepted for analysis.
%       aux  - <struct> .rippMaps .rippStates (accepted-based rate table) .sig
%                       (signal bundle, reuse across methods) .spkTimes .muTimes
%                       .uType .boutTimes .vldTimes .nremTimes .fs .rippCh.
%
%   DEPENDENCIES:
%       basepaths2vars, ripp_methods, evt_spkPrep, evt_boutTimes, ripp_pickCh,
%       ripp_sigLoad, ripp_sigPrep, ripp_noiseFloor, ripp_times, ripp_params,
%       evt_maps, evt_emgScore, evt_spkGain, evt_ctrlTimes, evt_states, evt_qa.
%
%   HISTORY:
%       260719 factored out of ripp_wrapper as the shared detection core.

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
addParameter(p, 'mapDur', [-0.1 0.1], @isnumeric);
addParameter(p, 'verbose', false, @islogical);
parse(p, basepath, varargin{:});
met     = p.Results.met;
win     = p.Results.win;
v       = p.Results.v;
rippCh  = p.Results.rippCh;
sig     = p.Results.sig;
mapDur  = p.Results.mapDur;
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
rippMaps = evt_maps(rippSig, ripp.peakTime, fs, 'mapDur', mapDur);

%% ========================================================================
%  QA METRICS + ACCEPTANCE (mark, do not remove)
%  ========================================================================

ripp.emg       = evt_emgScore(emg, ripp.times, fs, 'baselineTimes', nremTimes);
ripp.spkGain   = evt_spkGain(muTimes, ripp.times);
ripp.ctrlTimes = evt_ctrlTimes(ripp.times, 'vldTimes', vldTimes, 'flgPlot', false);

% accepted = in-state (NREM or valid) AND low EMG AND above the MUA-gain gate
if met.nremOnly, inTimes = nremTimes; else, inTimes = vldTimes; end
ripp.accepted = evt_qa(ripp.peakTime, 'inTimes', inTimes, ...
    'metrics', [ripp.emg, ripp.spkGain], ...
    'ranges', {[-Inf, met.thrEmg], [met.gainThr, Inf]}, ...
    'names', {'emg', 'gain'});

% vigilance state per event (all events); rate/density over accepted only
[ripp.state, rippStates] = evt_states(ripp.times, ripp.peakTime, boutTimes, ...
    'accepted', ripp.accepted, 'basepath', basepath, ...
    'flgSave', false, 'flgPlot', false, 'name', 'ripp', 'lbl', 'Ripple');

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

aux.rippMaps   = rippMaps;
aux.rippStates = rippStates;
aux.sig        = sig;
aux.spkTimes   = spkTimes;
aux.muTimes    = muTimes;
aux.uType      = uType;
aux.boutTimes  = boutTimes;
aux.vldTimes   = vldTimes;
aux.nremTimes  = nremTimes;
aux.fs         = fs;
aux.rippCh     = rippCh;

end     % EOF
