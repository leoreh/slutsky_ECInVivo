function ch = ed_pickCh(basepath, varargin)
% ED_PICKCH Resolve the ED detection channel (1-indexed).
%
%   ch = ED_PICKCH(basepath, varargin)
%
%   SUMMARY:
%       The ED twin of ripp_pickCh, and it exists because the two events want
%       different electrodes: ripples live in the CA1 pyramidal layer and are
%       found by ripple-band POWER, whereas a discharge is a rare, solitary
%       transient. Two steps:
%
%         1. DROP a channel whose crossing rate is a Tukey outlier among the
%            probe's channels (above Q3 + 1.5 IQR of the fraction of samples
%            over the detection threshold). Such a channel is not a better
%            electrode, it is a broken one - in raMCU5 it produced 13457
%            candidates where every other channel gave 4760-6495.
%         2. PICK, among the rest, the largest transients in the band:
%            prctile(|filt|, 99.9). A discharge is a large transient, so this
%            is the thing itself rather than a proxy for it.
%
%       Both steps were measured against ground truth - the detector run on
%       EVERY channel of raMCU3/4/5, scored by how well fastZ separates the
%       curated discharges (dev/ed_chSweep.m). The rule lands on AUC 0.957 /
%       0.986 / 0.978 where the best channel gives 0.960 / 0.986 / 0.983.
%       Without step 1 raMCU5 falls to 0.951 and doubles the curation load.
%
%       Two alternatives were tried and refuted. Normalising the percentile by
%       the channel's own background rewards whichever channel is quietest.
%       Ranking by FEWEST candidates - the intuition that every channel sees
%       the discharge so the cleanest wins - is anti-correlated with AUC in
%       raMCU3 (-0.69): a channel with few crossings is one where the
%       discharge itself barely crosses.
%
%       Resolves in priority:
%         1. an explicit 'edCh' (caller override);
%         2. the channel the pipeline already detected on, read from
%            <basename>.ed.mat (ed.info.edCh), so a re-run and every consumer
%            stay on one electrode;
%         3. the producer pick - the score above, over probe windows spread
%            across the recording.
%
%       ONE channel, never an average. Averaging is what the sleep_sig eeg did,
%       and it mixed sites that differ: in raMCU1 it spanned two shanks, in
%       raMCU2 it included the weakest channel of fifteen. A discharge is a
%       laminar event, so averaging across sites attenuates it by an amount
%       that varies per mouse - which is not a difference in the biology.
%
%       Only channels that belong to a spike group are considered; the rest are
%       auxiliary (accelerometer and the like) and carry no LFP.
%
%   INPUTS:
%       basepath - <char> session directory.
%       varargin - Parameter/Value:
%           'basename' - <char>   file stem. {folder name}
%           'session'  - <struct> session metadata. {loaded if needed}
%           'edCh'     - <num>    explicit channel; overrides all else. {[]}
%           'win'      - <vec>    window [start end] (s) to probe. {[0 Inf]}
%           'passband' - <vec>    band the score is measured in (Hz). {[60 150]}
%           'thr'      - <num>    threshold the crossing rate uses. {8}
%           'flgForce' - <log>    skip ed.mat; force the producer pick.{false}
%
%   OUTPUT:
%       ch       - <num>  1-indexed detection channel.
%
%   DEPENDENCIES:
%       binary_load, filterLFP, basepaths2vars.
%
%   HISTORY:
%       260721 created. Replaces the sleep_sig eeg average, which was a
%              different set of channels in every mouse (see
%              dev/ed_pipeline_rebuild.md).

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'basepath', @ischar);
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'session', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'edCh', [], @isnumeric);
addParameter(p, 'win', [0 Inf], @isnumeric);
addParameter(p, 'passband', [60 150], @isnumeric);
addParameter(p, 'thr', 8, @isnumeric);
addParameter(p, 'flgForce', false, @islogical);
parse(p, basepath, varargin{:});
session  = p.Results.session;
edCh     = p.Results.edCh;
win      = p.Results.win;
passband = p.Results.passband;
thr      = p.Results.thr;
flgForce = p.Results.flgForce;

basename = p.Results.basename;
if isempty(basename), [~, basename] = fileparts(basepath); end

% 1. explicit override
if ~isempty(edCh)
    ch = edCh;
    return;
end

% 2. consumer: the channel detection already ran on
if ~flgForce
    f = fullfile(basepath, [basename, '.ed.mat']);
    if isfile(f)
        S = load(f, 'ed');
        if isfield(S, 'ed') && isfield(S.ed, 'info') ...
                && isfield(S.ed.info, 'edCh') && ~isempty(S.ed.info.edCh)
            ch = S.ed.info.edCh;
            return;
        end
    end
end

% 3. producer pick
if isempty(session)
    v = basepaths2vars('basepaths', {basepath}, 'vars', {'session'}, ...
        'flgPrnt', false);
    session = v.session;
end
ch = bestEdCh(basepath, basename, session, win, passband, thr);

end     % EOF


% =========================================================================
%  LOCAL
% =========================================================================
function ch = bestEdCh(basepath, basename, session, win, passband, thr)
% Largest transients in the detection band, among channels that are not
% pathologically noisy. Probe windows are SPREAD across the recording rather
% than contiguous, because a single block can land entirely in one vigilance
% state and a discharge is a sleep event.
NPROBE  = 20;               % probe windows
PROBDUR = 15;               % duration of each [s]
PRC     = 99.9;             % percentile standing for "the biggest transients"

fs   = session.extracellular.srLfp;
nCh  = session.extracellular.nChannels;
chOk = sort(unique([session.extracellular.spikeGroups.channels{:}]));
if round(session.extracellular.sr) == 24414, b2u = 1; else, b2u = 0.195; end

fname = fullfile(basepath, [basename, '.lfp']);
d = dir(fname);
if isempty(d)
    error('ed_pickCh:noLfp', '%s not found', fname);
end
durRec = d.bytes / 2 / nCh / fs;

t1 = max(0, win(1));
if isinf(win(2)), t2 = durRec; else, t2 = min(win(2), durRec); end
tProbe = linspace(t1, t2 - PROBDUR, NPROBE);

sig = cell(NPROBE, 1);
for iPrb = 1 : NPROBE
    sig{iPrb} = double(binary_load(fname, 'fs', fs, 'nCh', nCh, ...
        'start', tProbe(iPrb), 'duration', PROBDUR, 'ch', chOk, ...
        'bit2uv', b2u));
end
sig = vertcat(sig{:});

score = zeros(1, numel(chOk));
xRate = zeros(1, numel(chOk));
for iCh = 1 : numel(chOk)
    fb = filterLFP(sig(:, iCh), 'fs', fs, 'type', 'butter', ...
        'dataOnly', true, 'order', 5, 'passband', passband, ...
        'graphics', false);
    scl = 1.4826 * median(abs(fb - median(fb)));
    score(iCh) = prctile(abs(fb), PRC);
    xRate(iCh) = mean(abs(fb) > thr * scl);
end

% step 1: drop the pathologically noisy channels (Tukey upper fence)
q = prctile(xRate, [25 75]);
ok = xRate <= q(2) + 1.5 * (q(2) - q(1));
if ~any(ok), ok = true(size(ok)); end

% step 2: largest transients among the rest
score(~ok) = -Inf;
[~, iBest] = max(score);
ch = chOk(iBest);

end     % bestEdCh
