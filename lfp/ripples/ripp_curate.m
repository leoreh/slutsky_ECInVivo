function [ripp, hFig] = ripp_curate(basepath, varargin)
% RIPP_CURATE Curate ripples by waveform TYPE (stage 2).
%
%   [ripp, hFig] = RIPP_CURATE(basepath, varargin)
%
%   SUMMARY:
%       The curation stage of the ripple pipeline (detect -> curate -> analyze).
%       It loads <basename>.ripp.mat and its per-event waveforms, applies
%       met.qa as a POOL, groups the pool into waveform clusters, and lets a
%       human accept or reject whole shapes. The GUI itself is evt_curate,
%       shared with the ED pipeline; everything here is loading and the two
%       ways to drive it:
%
%       - Headless (flgGui = false): apply met.qa, save .accepted (and the spec
%         in ripp.info.qa) back to ripp.mat, rebuild rippStates and clear the
%         stale analyze products. This is the automatic gate - a batch run needs
%         no human, and no waveforms are loaded.
%
%       - Interactive (flgGui = true, default): the cluster GUI. Read
%         evt_curate for what the controls do.
%
%       WHY SHAPES AND NOT THRESHOLDS. The metric gate cannot see what an event
%       looks like. It removes what is loud in the EMG and what nothing fires
%       during, and that is all it can do - a step artifact with quiet muscle
%       and a bystander burst passes every threshold there is. Blind waveform
%       clustering separates a population of steps from a population of ripples
%       in one decision, which is what a session with movement artifact needs,
%       and it costs a dozen judgements instead of twenty thousand.
%
%       WHAT IS CLUSTERED is the raw LFP around the peak; the band-passed trace
%       rides along as a second Y option, because whether a cluster actually
%       oscillates is usually the thing that settles it. A ripple separates from
%       a transient inside +-30 ms (met.clust.win), which is why the window is
%       narrower than the ED pipeline's +-50 ms - there, the DECAY is the tell.
%
%   INPUTS:
%       basepath - <char> session directory (must hold <basename>.ripp.mat).
%       varargin - Parameter/Value:
%           'basename' - <char>   file stem. {folder name}
%           'met'      - <struct> config; reads .qa and .clust.
%                                 {ripp_methods('default')}
%           'flgGui'   - <log>    open the GUI (true) or gate headless. {true}
%           'flgInvalidate' - <log> when the mask changes, delete the stale
%                                 accepted-aligned analyze products (rippSpks /
%                                 rippSpkMaps / rippSpkLfp) so they cannot be
%                                 read stale before ripp_analyze reruns. The
%                                 batch turns this off (analyze overwrites them
%                                 next). {true}
%           'Visible'  - <char>   'on' | 'off' for headless GUI tests. {'on'}
%           'verbose'  - <log>    print progress? {true}
%
%   OUTPUTS:
%       ripp - <struct> the loaded events, with .accepted as of the call.
%       hFig - <handle> the GUI figure ([] when headless).
%
%   DEPENDENCIES:
%       evt_files, evt_curate, evt_detrend, ripp_methods, ripp_invalidate;
%       maps rebuild: basepaths2vars, evt_boutTimes, ripp_sigLoad,
%       ripp_sigPrep, evt_maps.
%
%   HISTORY:
%       260719b the curation stage; absorbed evt_qa's ripple role.
%       260720  waveform view tiled by state; maps READ from rippMaps.
%       260722  rebuilt on waveform clustering, sharing evt_curate with the ED
%               pipeline - accept TYPES, not threshold values. The threshold
%               knobs are now generated from met.qa.ranges, so the gate is
%               described once (ripp_methods) instead of twice. 'qa' became
%               'met', because clustering needs met.clust.

%% ========================================================================
%  ARGUMENTS + LOAD
%  ========================================================================
p = inputParser;
addRequired(p, 'basepath', @ischar);
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'met', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'flgGui', true, @islogical);
addParameter(p, 'flgInvalidate', true, @islogical);
addParameter(p, 'Visible', 'on', @(x) any(strcmpi(char(x), {'on', 'off'})));
addParameter(p, 'verbose', true, @islogical);
parse(p, basepath, varargin{:});
met      = p.Results.met;
flgGui   = p.Results.flgGui;
flgInval = p.Results.flgInvalidate;
verbose  = p.Results.verbose;

basename = p.Results.basename;
if isempty(basename), [~, basename] = fileparts(basepath); end
if isempty(met), met = ripp_methods('default'); end

files = evt_files(basepath, basename, 'ripp');
if ~isfile(files.evt)
    error('ripp_curate:noFile', ...
        '%s not found; run detection first (ripp_wrapper).', files.evt);
end
S = load(files.evt, 'ripp');
ripp = S.ripp;

cfg = struct('met', met, 'file', files.evt, 'var', 'ripp', ...
    'basepath', basepath, 'basename', basename, 'lbl', 'Ripple', ...
    'flgGui', flgGui, 'Visible', char(p.Results.Visible), ...
    'onSaved', @(changed) invalidate(changed, flgInval, basepath, basename, ...
    verbose));

%% ========================================================================
%  CURATE
%  ========================================================================
% Waveforms are the clustering's whole input, so they are loaded only when
% there is a human to look at them - the headless gate is a batch operation and
% must not pay for a 50 MB read it has no use for.
wv = struct();
tst = [];
if flgGui
    [wv, tst] = loadMaps(basepath, basename, ripp, files.maps);
end

[ripp.accepted, hFig] = evt_curate(ripp, wv, tst, cfg);

% Only headless has a final answer to report. In the GUI this mask is the seed
% the window opened on, and printing it as a result would be a lie the moment
% the user ticks anything.
if verbose && ~flgGui
    fprintf('[RIPP_CURATE] %s : %d / %d accepted (headless)\n', basename, ...
        nnz(ripp.accepted), numel(ripp.accepted));
end

end     % EOF


% =========================================================================
%  LOCAL
% =========================================================================
function msg = invalidate(changed, flgInval, basepath, basename, verbose)
% Drop the analyze products a moved mask has staled. Returns what to say about
% it, which evt_curate appends to the GUI's save notification ('' = nothing to
% say); printed as well, because the headless path has no notification and
% ripp_invalidate itself deletes silently.
msg = '';
if ~changed || ~flgInval, return; end
nDel = ripp_invalidate(basepath, basename);
if nDel > 0
    msg = sprintf('Removed %d stale product(s) - rerun ripp_analyze.', nDel);
    if verbose, fprintf('[RIPP_CURATE] %s\n', msg); end
end

end     % invalidate


function [wv, tst] = loadMaps(basepath, basename, ripp, fileMaps)
% Per-event waveforms behind the clustering and the view: the raw LFP and the
% band-passed trace, detrended, cropped to DISPDUR.
%
% DETREND BEFORE CROP. evt_detrend fits its baseline on the flanks of whatever
% window it is handed, and on a +-60 ms crop the flanks are still inside the
% sharp wave - so the line would be fitted on the event itself. The saved map is
% +-100 ms, wide enough for the flanks to be baseline.
%
% Detection writes rippMaps over ALL detected events (ripp_wrapper), so this is
% normally a file read. A session detected before that convention, or one whose
% file no longer matches the event list, falls back to rebuilding the detection
% signal exactly as detection did.
dispDur = [-0.06 0.06];
FLDS = {'lfp', 'filt'};             % clustered on the first, both viewable

maps = readMaps(fileMaps, numel(ripp.peakTime));
if isempty(maps)
    maps = rebuildMaps(basepath, basename, ripp, dispDur);
end

tst = maps.tstamps;
keep = tst >= dispDur(1) & tst <= dispDur(2);
wv = struct();
for iFld = 1 : numel(FLDS)
    if ~isfield(maps, FLDS{iFld}), continue; end
    % kept SINGLE, as the file stores them: 100k events x 151 samples x two
    % traces is 137 MB single and 274 MB double, and nothing downstream needs
    % the precision - evt_clust casts its own window and guiTbl_xy plots either
    m = evt_detrend(double(maps.(FLDS{iFld})), tst);
    wv.(FLDS{iFld}) = single(m(:, keep));
    clear m                 % 215 MB, before the next field allocates its own
end
tst = tst(keep);
if isempty(fieldnames(wv))
    error('ripp_curate:noMaps', ...
        'rippMaps carries no lfp; re-run detection with flgSave.');
end

end     % loadMaps


function maps = readMaps(file, nEv)
% The saved maps, or [] when they are absent or belong to other events.
maps = [];
if ~isfile(file), return; end
S = load(file, 'rippMaps');
if ~isfield(S, 'rippMaps') || ~isfield(S.rippMaps, 'lfp'), return; end
if size(S.rippMaps.lfp, 1) ~= nEv, return; end
maps = S.rippMaps;

end     % readMaps


function maps = rebuildMaps(basepath, basename, ripp, dispDur)
% Rebuild the detection signal exactly as detection did, then cut the maps from
% it. The window is the recording frame ripp.info.win, so the absolute event
% times line up with the signal. Cut twice DISPDUR so evt_detrend still has
% flanks outside the event to fit on.
win = [0 Inf];
if isfield(ripp, 'info') && isfield(ripp.info, 'win'), win = ripp.info.win; end

v = basepaths2vars('basepaths', {basepath}, ...
    'vars', {'session', 'sleep_states'});
fs = v.session.extracellular.srLfp;
if isinf(win(2)), win(2) = v.session.extracellular.nSamples / fs; end
w0 = win(1);
if ~isfinite(w0), w0 = 0; end

[~, ~, nremTimes] = evt_boutTimes(v, win, win(2) - win(1));
lfp = ripp_sigLoad(basepath, 'win', win, 'session', v.session, ...
    'basename', basename, 'rippCh', ripp.info.rippCh, 'bit2uv', []);
% a pre-260720 ripp.mat carries no otlThr, hence the default
otlThr = 8;
if isfield(ripp.info, 'otlThr'), otlThr = ripp.info.otlThr; end
rippSig = ripp_sigPrep(lfp, fs, 'detectMet', ripp.info.detectMet, ...
    'passband', ripp.info.passband, 'zMet', ripp.info.zMet, ...
    'nremTimes', nremTimes, 'otlThr', otlThr);

maps = evt_maps(rippSig, ripp.peakTime - w0, fs, 'mapDur', 2 * dispDur);

end     % rebuildMaps
