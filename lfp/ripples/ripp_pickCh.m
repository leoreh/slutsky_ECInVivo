function ch = ripp_pickCh(basepath, varargin)
% RIPP_PICKCH Resolve the ripple detection channel (1-indexed).
%
%   ch = RIPP_PICKCH(basepath, varargin)
%
%   SUMMARY:
%       The single ripple-owned channel resolver. Replaces evt_pickCh and
%       evt_rippCh and absorbs the screen's best-channel probe, retiring
%       session.channelTags.Ripple. Resolves in priority:
%         1. an explicit 'rippCh' (caller override);
%         2. the channel the pipeline already detected on, read from
%            <basename>.ripp.mat (ripp.info.rippCh) - the single source of truth,
%            so consumers (ED, the curation GUI) follow the ripple analysis;
%         3. the producer pick - the channel with the most ripple-band power in a
%            short NREM probe, when a session is available and no ripple output
%            exists yet.
%       Falls back to channel 1 with a warning when none of these can be formed
%       (a consumer asked before ripples were detected and without a session).
%       Channels are 1-indexed throughout (the binary_load convention).
%
%   INPUTS:
%       basepath - <char> session directory.
%       varargin - Parameter/Value:
%           'basename'  - <char>   file stem. {folder name}
%           'session'   - <struct> session metadata (for the producer pick). {[]}
%           'rippCh'    - <num>    explicit channel(s); overrides all else. {[]}
%           'win'       - <vec>    window [start end] (s) for the probe. {[0 Inf]}
%           'nremTimes' - <mat>    [N x 2] NREM bouts (s) to place the probe. {[]}
%           'flgForce'  - <log>    skip ripp.mat; force the producer pick. {false}
%           'passband'  - <vec>    ripple band for the probe (Hz). {[120 220]}
%
%   OUTPUT:
%       ch       - <num>    1-indexed detection channel(s).
%
%   DEPENDENCIES:
%       binary_load, filterLFP; basepaths2vars (only to load a session for the
%       producer pick when one is not supplied).
%
%   HISTORY:
%       260719 merge evt_pickCh + evt_rippCh + the screen's bestRippCh; retire
%              channelTags.Ripple.

p = inputParser;
addRequired(p, 'basepath', @ischar);
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'session', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'rippCh', [], @isnumeric);
addParameter(p, 'win', [0 Inf], @isnumeric);
addParameter(p, 'nremTimes', [], @isnumeric);
addParameter(p, 'flgForce', false, @islogical);
addParameter(p, 'passband', [120 220], @isnumeric);
parse(p, basepath, varargin{:});
basename  = p.Results.basename;
session   = p.Results.session;
rippCh    = p.Results.rippCh;
win       = p.Results.win;
nremTimes = p.Results.nremTimes;
flgForce  = p.Results.flgForce;
passband  = p.Results.passband;

if isempty(basename)
    [~, basename] = fileparts(basepath);
end

% 1. explicit override
if ~isempty(rippCh)
    ch = rippCh;
    return;
end

% 2. consumer: the channel ripples detected on (the single source of truth)
if ~flgForce
    f = fullfile(basepath, [basename, '.ripp.mat']);
    if isfile(f)
        S = load(f, 'ripp');
        if isfield(S, 'ripp') && isfield(S.ripp, 'info') ...
                && isfield(S.ripp.info, 'rippCh') && ~isempty(S.ripp.info.rippCh)
            ch = S.ripp.info.rippCh;
            return;
        end
    end
end

% 3. producer: best ripple-band channel in a short NREM probe
if isempty(session)
    v = basepaths2vars('basepaths', {basepath}, 'vars', {'session'}, ...
        'flgPrnt', false);
    if isfield(v, 'session'), session = v.session; end
end
hasExtra = ~isempty(session) && isfield(session, 'extracellular') ...
    && isfield(session.extracellular, 'srLfp') ...
    && isfield(session.extracellular, 'nChannels');
if hasExtra && isfile(fullfile(basepath, [basename, '.lfp']))
    ch = bestRippCh(basepath, basename, session, win, nremTimes, passband);
    return;
end

% 4. last resort
warning('ripp_pickCh:noChannel', ...
    'no ripple channel for %s; using channel 1.', basename);
ch = 1;

end     % EOF


% =========================================================================
%  LOCALS
% =========================================================================
function ch = bestRippCh(basepath, basename, session, win, nremTimes, passband)
% channel with the most ripple-band power in a short NREM probe (1-indexed)
fs  = session.extracellular.srLfp;
nCh = session.extracellular.nChannels;
if round(session.extracellular.sr) == 24414, b2u = 1; else, b2u = 0.195; end

% probe span: the first NREM bout > 30 s, else a window at the start
if isinf(win(2)), winEnd = win(1) + 120; else, winEnd = win(2); end
if ~isempty(nremTimes)
    idx = find(diff(nremTimes, [], 2) > 30, 1);
    if isempty(idx), idx = 1; end
    pT = nremTimes(idx, :) + win(1);
else
    pT = [win(1), win(1) + min(120, winEnd - win(1))];
end
dur = min(120, pT(2) - pT(1));

probe = double(binary_load(fullfile(basepath, [basename, '.lfp']), ...
    'fs', fs, 'nCh', nCh, 'start', pT(1), 'duration', dur, ...
    'ch', 1:nCh, 'bit2uv', b2u));
rmsCh = zeros(1, nCh);
for iCh = 1:nCh
    fb = filterLFP(probe(:, iCh), 'fs', fs, 'type', 'butter', ...
        'dataOnly', true, 'order', 5, 'passband', passband, 'graphics', false);
    rmsCh(iCh) = rms(fb);
end
[~, ch] = max(rmsCh);
end
