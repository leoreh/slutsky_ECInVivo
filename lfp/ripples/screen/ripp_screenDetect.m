function [ripp, rippStates, meta] = ripp_screenDetect(basepath, method, varargin)
% RIPP_SCREENDETECT Run one detection method on one session slice (light path).
%
%   [ripp, rippStates, meta] = RIPP_SCREENDETECT(basepath, method, varargin)
%
%   SUMMARY:
%       The detection half of the screen. Runs the existing ripple pipeline
%       pieces - signal load, prep, threshold, per-event params, state labels -
%       for a single method (from ripp_screenMethods) over a window, and stops
%       before the heavy spike/map/phase steps. Writes NOTHING (evt_states is
%       called with flgSave=false), so the canonical <basename>.ripp.mat is
%       never touched. Also scores each event against multi-unit spiking (MUA)
%       as a label-free detector-quality metric.
%
%   INPUTS:
%       basepath - (Char)   session directory.
%       method   - (Struct) one element of ripp_screenMethods().
%       varargin - Parameter/Value:
%           'win'     - (Vec)    window [start end] (s). {[0 Inf]}
%           'v'       - (Struct) pre-loaded basepaths2vars struct (session,
%                                sleep_states, spikes). Loaded if empty.
%           'verbose' - (Log)    print progress. {false}
%
%   OUTPUTS:
%       ripp       - (Struct) events + per-event params (.times .peakTime .amp
%                             .freq .freqEvent .freqPeak .energy .dur .skew
%                             .state .accepted .muaZ .info). Times ABSOLUTE.
%       rippStates - (Table)  per-bout rate/density (evt_states 2nd output).
%       meta       - (Struct) .rippCh .fs .win .nEvents .rateHz .muaPos .method.
%
%   DEPENDENCIES:
%       basepaths2vars, evt_boutTimes, ripp_sigLoad, ripp_sigPrep, ripp_times,
%       ripp_params, evt_states, binary_load, filterLFP.
%
%   HISTORY:
%       Created: 260716 (detection-review parameter screen).

%% ---- arguments ----
p = inputParser;
addRequired(p, 'basepath', @ischar);
addRequired(p, 'method', @isstruct);
addParameter(p, 'win', [0 Inf], @isnumeric);
addParameter(p, 'v', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'verbose', false, @islogical);
parse(p, basepath, method, varargin{:});
win     = p.Results.win;
v       = p.Results.v;
verbose = p.Results.verbose;
[~, basename] = fileparts(basepath);

%% ---- session + optional data ----
% basepaths2vars now matches each file by its dotted name, so 'sleep_states'
% no longer picks up AccuSleep_states.mat (see io/basepaths2vars.m).
if isempty(v)
    v = basepaths2vars('basepaths', {basepath}, ...
        'vars', {'session', 'sleep_states', 'spikes'});
end
ses = v.session;
fs  = ses.extracellular.srLfp;
nCh = ses.extracellular.nChannels;
if round(ses.extracellular.sr) == 24414, b2u = 1; else, b2u = 0.195; end
if isinf(win(2)), win(2) = ses.extracellular.nSamples / fs; end
sigDur = win(2) - win(1);

% window-relative bout / NREM times (empty degrades gracefully downstream)
[boutTimes, ~, nremTimes] = evt_boutTimes(v, win, sigDur);

%% ---- detection channel ----
switch method.chMode
    case 'tag'
        rippCh = evt_pickChSafe(ses);           % averaged by ripp_sigLoad
    case 'best'
        rippCh = bestRippCh(basepath, basename, win, nremTimes, fs, nCh, b2u);
    otherwise
        error('ripp_screenDetect:chMode', 'unknown chMode "%s"', method.chMode);
end

%% ---- light detection path ----
if verbose, fprintf('[SCREEN] %s / %s : ch %s\n', basename, ...
        method.name, mat2str(rippCh)); end
[lfp, ~, fs] = ripp_sigLoad(basepath, 'win', win, 'session', ses, ...
    'basename', basename, 'rippCh', rippCh, 'bit2uv', []);
rippSig = ripp_sigPrep(lfp, fs, 'detectMet', method.detectMet, ...
    'passband', method.passband, 'zMet', method.zMet, 'nremTimes', nremTimes);
ripp = ripp_times(rippSig, fs, 'thr', method.thr, 'limDur', method.limDur);
ripp = ripp_params(rippSig, ripp);
[ripp.state, rippStates] = evt_states(ripp.times, ripp.peakTime, boutTimes, ...
    'basepath', basepath, 'flgPlot', false, 'flgSave', false, ...
    'name', 'ripp', 'lbl', 'Ripple');
ripp.accepted = true(size(ripp.times, 1), 1);

%% ---- MUA convergence (label-free detector quality) ----
mua = [];
if isfield(v, 'spikes') && isfield(v.spikes, 'times')
    mua = sort(vertcat(v.spikes.times{:})) - win(1);
    mua = mua(mua > 0 & mua < sigDur);
end
ripp.muaZ = muaGain(mua, ripp.peakTime, nremTimes, sigDur);

%% ---- absolute time + provenance ----
ripp.times    = ripp.times + win(1);
ripp.peakTime = ripp.peakTime + win(1);
ripp.info.basename = basename;
ripp.info.rippCh   = rippCh;
ripp.info.passband = method.passband;
ripp.info.zMet     = method.zMet;
ripp.info.method   = method.name;
ripp.info.win      = win;

meta = struct('rippCh', rippCh, 'fs', fs, 'win', win, ...
    'nEvents', size(ripp.times, 1), 'rateHz', size(ripp.times, 1) / sigDur, ...
    'muaPos', mean(ripp.muaZ > 1), 'method', method.name);

end     % EOF


% =========================================================================
%  LOCALS
% =========================================================================
function ch = evt_pickChSafe(ses)
ch = 1;
if isfield(ses, 'channelTags') && isfield(ses.channelTags, 'Ripple') ...
        && ~isempty(ses.channelTags.Ripple)
    ch = ses.channelTags.Ripple;
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
for c = 1:nCh
    fb = filterLFP(probe(:, c), 'fs', fs, 'type', 'butter', 'dataOnly', true, ...
        'order', 5, 'passband', [120 220], 'graphics', false);
    rmsCh(c) = rms(fb);
end
[~, ch] = max(rmsCh);
end

% -------------------------------------------------------------------------
function z = muaGain(mua, peakTime, nremTimes, sigDur)
% per-event MUA z: spikes in +/-25 ms vs the NREM 50-ms-bin rate distribution
z = nan(numel(peakTime), 1);
if isempty(mua), return; end
half  = 0.025;
edges = (0 : 2*half : sigDur)';
cnt   = histcounts(mua, edges)';
bc    = edges(1:end-1) + half;

% baseline bins = those whose centre falls in NREM (all bins if no NREM)
if isempty(nremTimes)
    nremBin = true(numel(cnt), 1);
else
    nremBin = false(numel(cnt), 1);
    for i = 1:size(nremTimes, 1)
        nremBin(bc >= nremTimes(i, 1) & bc <= nremTimes(i, 2)) = true;
    end
    if ~any(nremBin), nremBin = true(numel(cnt), 1); end
end
mu = mean(cnt(nremBin), 'omitnan');
sd = std(cnt(nremBin), 'omitnan');
if sd == 0, sd = 1; end

for i = 1:numel(peakTime)
    n = sum(mua >= peakTime(i) - half & mua < peakTime(i) + half);
    z(i) = (n - mu) / sd;
end
end
