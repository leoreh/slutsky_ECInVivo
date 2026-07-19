function [det, meta] = ripp_screenDetect(basepath, methods, varargin)
% RIPP_SCREENDETECT Run every method on one session slice (load once).
%
%   [det, meta] = RIPP_SCREENDETECT(basepath, methods, varargin)
%
%   SUMMARY:
%       The detection half of the screen, for one session. Loads the signals
%       once, prepares the filtered detection signal once per unique signal
%       configuration (chMode, passband, detectMet, zMet), and then runs every
%       method on the shared signal. Methods that differ only in the threshold
%       therefore reuse the same filtered trace instead of reloading and
%       refiltering. Runs the existing pieces (ripp_sigLoad, ripp_sigPrep,
%       ripp_times, ripp_params, evt_states) and stops before the heavy
%       spike/map/phase steps. Writes NOTHING - evt_states is called with
%       flgSave=false, so <basename>.ripp.mat is never touched. Scores each
%       event against multi-unit spiking (MUA) as a label-free quality metric.
%       A method with .calibThr replaces its fixed peak threshold with one set
%       to the 1/f noise floor (ripp_noiseFloor).
%
%   INPUTS:
%       basepath - <char>   session directory.
%       methods  - <struct> ripp_screenMethods() array.
%       varargin - Parameter/Value:
%           'win'     - <vec>    window [start end] (s). {[0 3*3600]}
%           'v'       - <struct> pre-loaded basepaths2vars (session,
%                                sleep_states, spikes). Loaded if empty.
%           'verbose' - <log>    print progress. {false}
%
%   OUTPUTS:
%       det  - <struct> [1 x nMethod] with .ripp (events + per-event params,
%                       absolute times), .rippStates (evt_states bout table),
%                       .name.
%       meta - <struct> [1 x nMethod] with .name .rippCh .fs .win .nEvents
%                       .rateHz .muaPos .thrPk .chi (chi/thrPk NaN unless
%                       calibrated).
%
%   DEPENDENCIES:
%       basepaths2vars, evt_boutTimes, evt_rippCh, ripp_sigLoad, ripp_sigPrep,
%       ripp_times, ripp_params, ripp_noiseFloor, evt_states, binary_load,
%       filterLFP.
%
%   HISTORY:
%       260716 detection-review parameter screen.
%       260717 load once; cache prepared signal across methods; calibThr path.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'basepath', @ischar);
addRequired(p, 'methods', @isstruct);
addParameter(p, 'win', [0 3 * 3600], @isnumeric);
addParameter(p, 'v', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'verbose', false, @islogical);
parse(p, basepath, methods, varargin{:});
win     = p.Results.win;
v       = p.Results.v;
verbose = p.Results.verbose;
[~, basename] = fileparts(basepath);
nM = numel(methods);

%% ========================================================================
%  SESSION + CONTEXT (loaded once)
%  ========================================================================

if isempty(v)
    v = basepaths2vars('basepaths', {basepath}, ...
        'vars', {'session', 'sleep_states', 'spikes'});
end
ses = v.session;
fs  = ses.extracellular.srLfp;
nCh = ses.extracellular.nChannels;
if round(ses.extracellular.sr) == 24414
    b2u = 1;
else
    b2u = 0.195;
end
if isinf(win(2))
    win(2) = ses.extracellular.nSamples / fs;
end
sigDur = win(2) - win(1);

% window-relative bout / NREM times (empty degrades gracefully downstream)
[boutTimes, ~, nremTimes] = evt_boutTimes(v, win, sigDur);

% multi-unit spike train, windowed (empty if no spikes)
mua = [];
if isfield(v, 'spikes') && isfield(v.spikes, 'times')
    mua = sort(vertcat(v.spikes.times{:})) - win(1);
    mua = mua(mua > 0 & mua < sigDur);
end

%% ========================================================================
%  DETECT (shared signal per config, threshold per method)
%  ========================================================================

det  = struct('ripp', cell(1, nM), 'rippStates', cell(1, nM), 'name', cell(1, nM));
meta = struct('name', cell(1, nM), 'rippCh', cell(1, nM), 'fs', cell(1, nM), ...
    'win', cell(1, nM), 'nEvents', cell(1, nM), 'rateHz', cell(1, nM), ...
    'muaPos', cell(1, nM), 'thrPk', cell(1, nM), 'chi', cell(1, nM));

sigCache = struct();
for iMethod = 1:nM
    method = methods(iMethod);

    % prepare the filtered signal once per unique configuration
    key = cfgKey(method);
    if ~isfield(sigCache, key)
        rippCh = pickChannel(method, basepath, basename, ses, win, ...
            nremTimes, fs, nCh, b2u);
        lfp = ripp_sigLoad(basepath, 'win', win, 'session', ses, ...
            'basename', basename, 'rippCh', rippCh, 'bit2uv', []);
        c.sig = ripp_sigPrep(lfp, fs, 'detectMet', method.detectMet, ...
            'passband', method.passband, 'zMet', method.zMet, ...
            'nremTimes', nremTimes);
        c.rippCh = rippCh;
        sigCache.(key) = c;
    end
    c = sigCache.(key);

    if verbose
        fprintf('[SCREEN] %s / %s : ch %s\n', basename, method.name, ...
            mat2str(c.rippCh));
    end

    % threshold: fixed, or calibrated to the 1/f noise floor
    thr = method.thr;
    chi = NaN;
    if method.calibThr
        [thrPk, nf] = ripp_noiseFloor(c.sig.lfp, fs, method, ...
            'nremTimes', nremTimes, 'targetFP', method.targetFP);
        thr(2) = thrPk;
        thr(1) = max(0.5, thrPk - (method.thr(2) - method.thr(1)));
        chi = nf.chi;
    end

    % detect + per-event params + state labels + MUA convergence
    ripp = ripp_times(c.sig, fs, 'thr', thr, 'limDur', method.limDur);
    ripp = ripp_params(c.sig, ripp);
    [ripp.state, rippStates] = evt_states(ripp.times, ripp.peakTime, ...
        boutTimes, 'basepath', basepath, 'flgPlot', false, ...
        'flgSave', false, 'name', 'ripp', 'lbl', 'Ripple');
    ripp.accepted = true(size(ripp.times, 1), 1);
    ripp.muaZ = muaGain(mua, ripp.peakTime, nremTimes, sigDur);

    % absolute time + provenance
    ripp.times     = ripp.times + win(1);
    ripp.peakTime  = ripp.peakTime + win(1);
    ripp.info.basename = basename;
    ripp.info.rippCh   = c.rippCh;
    ripp.info.passband = method.passband;
    ripp.info.method   = method.name;
    ripp.info.thr      = thr;
    ripp.info.win      = win;

    det(iMethod).ripp       = ripp;
    det(iMethod).rippStates = rippStates;
    det(iMethod).name       = method.name;

    meta(iMethod).name    = method.name;
    meta(iMethod).rippCh  = c.rippCh;
    meta(iMethod).fs      = fs;
    meta(iMethod).win     = win;
    meta(iMethod).nEvents = size(ripp.times, 1);
    meta(iMethod).rateHz  = size(ripp.times, 1) / sigDur;
    meta(iMethod).muaPos  = mean(ripp.muaZ > 1);
    meta(iMethod).thrPk   = thr(2);
    meta(iMethod).chi     = chi;
end

end     % EOF


% =========================================================================
%  LOCALS
% =========================================================================
function key = cfgKey(method)
% signal-config identity: methods with the same key share one prepared signal
key = matlab.lang.makeValidName(sprintf('%s_%s_%d_%s', method.chMode, ...
    num2str(method.passband), method.detectMet, method.zMet));
end

% -------------------------------------------------------------------------
function rippCh = pickChannel(method, basepath, basename, ses, win, ...
    nremTimes, fs, nCh, b2u)
% detection channel per chMode: the ripple tag, or the best NREM channel
switch method.chMode
    case 'tag'
        rippCh = evt_rippCh(basepath, basename, ses);
    case 'best'
        rippCh = bestRippCh(basepath, basename, win, nremTimes, fs, nCh, b2u);
    otherwise
        error('ripp_screenDetect:chMode', 'unknown chMode "%s"', method.chMode);
end
end

% -------------------------------------------------------------------------
function ch = bestRippCh(basepath, basename, win, nremTimes, fs, nCh, b2u)
% single channel with the most 120-220 Hz power in a short NREM probe
lfpFile = fullfile(basepath, [basename '.lfp']);
if ~isempty(nremTimes)
    idx = find(diff(nremTimes, [], 2) > 30, 1);   % first NREM bout > 30 s
    if isempty(idx), idx = 1; end
    pT = nremTimes(idx, :) + win(1);
else
    pT = [win(1), win(1) + min(120, win(2) - win(1))];
end
dur = min(120, pT(2) - pT(1));
probe = double(binary_load(lfpFile, 'fs', fs, 'nCh', nCh, ...
    'start', pT(1), 'duration', dur, 'ch', 1:nCh, 'bit2uv', b2u));
rmsCh = zeros(1, nCh);
for iCh = 1:nCh
    fb = filterLFP(probe(:, iCh), 'fs', fs, 'type', 'butter', ...
        'dataOnly', true, 'order', 5, 'passband', [120 220], ...
        'graphics', false);
    rmsCh(iCh) = rms(fb);
end
[~, ch] = max(rmsCh);
end

% -------------------------------------------------------------------------
function z = muaGain(mua, peakTime, nremTimes, sigDur)
% per-event MUA z: spikes in +/-25 ms vs the NREM 50-ms-bin rate distribution
z = nan(numel(peakTime), 1);
if isempty(mua), return; end
half  = 0.025;
edges = (0 : 2 * half : sigDur)';
cnt   = histcounts(mua, edges)';
bc    = edges(1:end-1) + half;

% baseline bins = those whose centre falls in NREM (all bins if no NREM)
if isempty(nremTimes)
    nremBin = true(numel(cnt), 1);
else
    nremBin = false(numel(cnt), 1);
    for iBout = 1:size(nremTimes, 1)
        nremBin(bc >= nremTimes(iBout, 1) & bc <= nremTimes(iBout, 2)) = true;
    end
    if ~any(nremBin), nremBin = true(numel(cnt), 1); end
end
mu = mean(cnt(nremBin), 'omitnan');
sd = std(cnt(nremBin), 'omitnan');
if sd == 0, sd = 1; end

for iEvt = 1:numel(peakTime)
    n = sum(mua >= peakTime(iEvt) - half & mua < peakTime(iEvt) + half);
    z(iEvt) = (n - mu) / sd;
end
end
