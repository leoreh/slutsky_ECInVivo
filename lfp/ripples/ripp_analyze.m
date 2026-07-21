function ripp = ripp_analyze(basepath, varargin)
% RIPP_ANALYZE Spike / phase analysis on the ACCEPTED ripples (stage 3).
%
%   ripp = RIPP_ANALYZE(basepath, varargin)
%
%   SUMMARY:
%       The heavy stage of the ripple pipeline (detect -> curate -> analyze). It
%       runs the expensive analyses only on the curated (accepted) events, so no
%       compute is wasted on the many events a per-mouse gate rejects. Loads the
%       saved <basename>.ripp.mat, takes the accepted subset, and computes: the
%       full spike modulation (evt_spks: per-unit stats, per-event population
%       metrics, 3D rasters, PETH) and spike-LFP phase coupling. Writes rippSpks
%       / rippSpkMaps / rippSpkLfp and folds the per-event population metrics
%       into ripp.spks. The per-event LFP maps are NOT built here - ripp_wrapper
%       writes them at detect over all events (see rippMaps there).
%
%       Run it after curation. In the batch (ripp_wrapper) the detect-stage
%       signal + spikes are passed through 'aux' so nothing is reloaded; run
%       standalone (after manual curation in ripp_curate), it reloads the signal
%       and spikes from disk. Either way it operates in the recording's frame via
%       ripp.info.win, so the saved absolute-time events line up with the signal.
%
%   INPUTS:
%       basepath - <char> session directory.
%       varargin - Parameter/Value:
%           'basename' - <char> file stem. {folder name}
%           'aux'      - <struct> detect-stage bundle (.sig .spkTimes .muTimes
%                                 .uType .fs) to skip the reload. {reload}
%           'mapDur'   - <vec>   PETH window [pre post] (s). {[-0.1 0.1]}
%           'flgSave'  - <log>   save the output files? {true}
%           'flgPlot'  - <log>   draw the spike-modulation summary? {true}
%           'verbose'  - <log>   print progress? {true}
%
%   OUTPUT:
%       ripp - <struct> the loaded events with .spks added (accepted-aligned).
%
%   FILES SAVED (when flgSave = true):
%       basename.ripp.mat        - re-saved with .spks (accepted-aligned)
%       basename.rippSpks.mat    - per-unit spike stats + PETH (light)
%       basename.rippSpkMaps.mat - 3D spike raster [unit x event x bin] (heavy)
%       basename.rippSpkLfp.mat  - spike-LFP phase coupling
%
%   DEPENDENCIES:
%       evt_files, evt_spks, evt_saveSpks, evt_plotSpks, spklfp_phase,
%       backup_file; reload path: basepaths2vars, evt_spkPrep, evt_boutTimes,
%       ripp_sigLoad, ripp_sigPrep.
%
%   HISTORY:
%       260719b split out of ripp_wrapper as the post-curation heavy stage.
%       260720  LFP maps moved to the detect stage (all events, mask-independent);
%               mapDur here now sizes the PETH only.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'basepath', @ischar);
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'aux', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'mapDur', [-0.1 0.1], @isnumeric);
addParameter(p, 'flgSave', true, @islogical);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'verbose', true, @islogical);
parse(p, basepath, varargin{:});
aux     = p.Results.aux;
mapDur  = p.Results.mapDur;
flgSave = p.Results.flgSave;
flgPlot = p.Results.flgPlot;
verbose = p.Results.verbose;

basename = p.Results.basename;
if isempty(basename), [~, basename] = fileparts(basepath); end
files = evt_files(basepath, basename, 'ripp');
files.phase = fullfile(basepath, [basename, '.rippSpkLfp.mat']);

S = load(files.evt, 'ripp');
ripp = S.ripp;
acc = ripp.accepted;

if ~any(acc)
    if verbose, fprintf('[RIPP_ANALYZE]: no accepted events; nothing to do.\n'); end
    return;
end

% recording frame: events are absolute, the signal starts at win(1)
win = [0 Inf];
if isfield(ripp, 'info') && isfield(ripp.info, 'win'), win = ripp.info.win; end
w0 = win(1);
if ~isfinite(w0), w0 = 0; end

%% ========================================================================
%  INPUTS (reuse the detect-stage bundle, or reload from disk)
%  ========================================================================

if ~isempty(aux)
    spkTimes = aux.spkTimes;
    muTimes  = aux.muTimes;
    uType    = aux.uType;
    rippSig  = aux.sig.rippSig;
    fs       = aux.fs;
else
    if verbose, fprintf('[RIPP_ANALYZE]: reloading signal + spikes...\n'); end
    v = basepaths2vars('basepaths', {basepath}, ...
        'vars', {'session', 'spikes', 'spktimes', 'sleep_states', 'units'});
    fs = v.session.extracellular.srLfp;
    if isinf(win(2)), win(2) = v.session.extracellular.nSamples / fs; end
    sigDur = win(2) - win(1);

    fsSpk = v.session.extracellular.sr;
    [spkTimes, muTimes, uType] = evt_spkPrep(v, win, sigDur, fsSpk);
    [~, ~, nremTimes] = evt_boutTimes(v, win, sigDur);

    lfp = ripp_sigLoad(basepath, 'win', win, 'session', v.session, ...
        'basename', basename, 'rippCh', ripp.info.rippCh, 'bit2uv', []);
    % rebuild the detection signal exactly as detection did, artifact mask
    % included; a pre-260720 ripp.mat carries no otlThr, hence the default
    otlThr = 8;
    if isfield(ripp.info, 'otlThr'), otlThr = ripp.info.otlThr; end
    rippSig = ripp_sigPrep(lfp, fs, 'detectMet', ripp.info.detectMet, ...
        'passband', ripp.info.passband, 'zMet', ripp.info.zMet, ...
        'nremTimes', nremTimes, 'otlThr', otlThr);
end

% window-relative times for the accepted events (signal starts at w0)
relTimes = ripp.times(acc, :)     - w0;
relCtrl  = ripp.ctrlTimes(acc, :) - w0;
relPeak  = ripp.peakTime(acc)     - w0;

%% ========================================================================
%  SPIKES + PHASE + MAPS (accepted events only)
%  ========================================================================

if verbose, fprintf('[RIPP_ANALYZE]: %d accepted events, analysing spikes...\n', nnz(acc)); end

rippSpks = evt_spks(spkTimes, muTimes, relTimes, relCtrl, relPeak, ...
    'unitType', uType, 'mapDur', mapDur);
ripp.spks = rippSpks.events;            % per-event population metrics (accepted)
rippSpks = rmfield(rippSpks, 'events');

spkLfp = spklfp_phase(rippSig.filt, spkTimes, fs, ...
    'lfpTimes', relTimes, 'nPerms', 0);

% the per-event LFP maps are NOT built here: ripp_wrapper writes them at detect
% over all events, row-aligned to ripp, and readers subset them by .accepted

%% ========================================================================
%  SAVE + PLOT
%  ========================================================================

if flgSave
    if verbose, fprintf('[RIPP_ANALYZE]: saving...\n'); end
    backup_file(files.evt);                      % preserve the curated mask
    save(files.evt, 'ripp', '-v7.3');           % re-save with .spks
    evt_saveSpks(rippSpks, basepath, basename, 'ripp');
    save(files.phase, 'spkLfp', '-v7.3');
end

if flgPlot
    evt_plotSpks(rippSpks, 'basepath', basepath, 'flgSaveFig', true, ...
        'name', 'ripp', 'lbl', 'Ripple');
end

if verbose, fprintf('[RIPP_ANALYZE]: done (%s).\n', basename); end

end     % EOF
