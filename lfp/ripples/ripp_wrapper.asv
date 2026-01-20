function ripp = ripp_wrapper(varargin)
% RIPP_WRAPPER Master function to run the complete ripple analysis pipeline.
%
% SUMMARY:
%   1. ripp_detect -> basename.ripp.mat
%   2. ripp_states -> basename.rippStates.mat
%   3. ripp_spks   -> basename.rippSpks.mat (with optional State restriction)
%   4. spklfp_calc -> (part of ripp struct in memory, optional save?)
%
% OUTPUT:
%   ripp - Combined structure for convenience (in memory).

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'rippCh', [], @isnumeric);
addParameter(p, 'thr', [1, 2.5, 1.2, 200, 100], @isnumeric);
addParameter(p, 'limDur', [20, 300, 20, 10], @isnumeric);
addParameter(p, 'passband', [100 300], @isnumeric);
addParameter(p, 'detectAlt', 3, @isnumeric);
addParameter(p, 'limState', [], @isnumeric);
addParameter(p, 'flgRefine', true, @islogical);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgSave', false, @islogical);
addParameter(p, 'flgSaveFig', false, @islogical);
addParameter(p, 'bit2uv', [], @isnumeric);
parse(p, varargin{:});

basepath = p.Results.basepath;
rippCh = p.Results.rippCh;
thr = p.Results.thr;
limDur = p.Results.limDur;
detectAlt = p.Results.detectAlt;
flgPlot = p.Results.flgPlot;
flgSave = p.Results.flgSave;
flgSaveFig = p.Results.flgSaveFig;
flgRefine = p.Results.flgRefine;
bit2uv = p.Results.bit2uv;
limState = p.Results.limState;

%% ========================================================================
%  SETUP
%  ========================================================================

cd(basepath);
[~, basename] = fileparts(basepath);
fprintf('--- Ripple Analysis for %s ---\n', basename);

% Load vars
vars = {'session', 'spikes', 'spktimes', 'sleep_states'};
v = basepaths2vars('basepaths', {basepath}, 'vars', vars);

% Load Session
fs = v.session.extracellular.srLfp;
fsSpk = v.session.extracellular.sr;
nchans = v.session.extracellular.nChannels;

% Ripple channel in LFP
if isempty(rippCh)
    if isfield(v.session.channelTags, 'Ripple')
        rippCh = v.session.channelTags.Ripple;
    else
        rippCh = 1;
    end
end

% Handle bit2uv
if isempty(bit2uv)
    if round(fsSpk) == 24414
        bit2uv = 1;
    else
        bit2uv = 0.195;
    end
end

% Load EMG
load([basename, '.sleep_sig.mat'], 'emg');

% Load LFP
fname = fullfile(basepath, [basename, '.lfp']);
sig = binary_load(fname, 'duration', Inf, 'fs', fs, 'nCh', nchans,...
    'start', 0, 'ch', rippCh, 'downsample', 1, 'bit2uv', bit2uv);

if size(sig, 2) > 1
    sig = mean(sig, 2);
end


%% ========================================================================
%  DETECT
%  ========================================================================
% Refinement run
if flgRefine
    winIdx = 1 : min(length(sig), 3 * 60 * 60 * fs);
    thrTmp = [1, 2.5, 1.2, 200, 50];

    rippTmp = ripp_detect('basepath', basepath, 'sig', sig(winIdx), 'emg', emg(winIdx),...
        'fs', fs, 'rippCh', rippCh, 'recWin', [0 Inf], 'thr', thrTmp, ...
        'detectAlt', detectAlt, 'flgPlot', false, 'flgSave', false);

    thr = round(rippTmp.info.thrData * 0.8, 1);
    thr(5) = 50;
    limDur = rippTmp.info.limDurData;
end

% Detection
ripp = ripp_detect('basepath', basepath, ...
    'sig', sig, ...
    'emg', emg, ...
    'fs', fs, ...
    'rippCh', rippCh, ...
    'recWin', [0 Inf], ...
    'thr', thr, ...
    'limDur', limDur, ...
    'detectAlt', detectAlt, ...
    'flgPlot', flgPlot, ...
    'flgSave', flgSave, ...
    'bit2uv', bit2uv);


%% ========================================================================
%  STATES
%  ========================================================================

try
    boutTimes = v.ss.bouts.times;
catch
    boutTimes = [];
end

ripp = ripp_states(ripp, boutTimes, ...
    'basepath', basepath, ...
    'flgPlot', flgPlot, ...
    'flgSave', flgSave, ...
    'flgSaveFig', flgSaveFig);


% Filter Ripples by State if requested
rippTimes = ripp.times;
if ~isempty(limState)
    stateIdx = ripp.states.idx(:, limState);
    rippTimes = ripp.times(stateIdx, :);
    fprintf('Limiting analysis to State %d\n', limState);
end


%% ========================================================================
%  SPIKES
%  ========================================================================

% Prep SU & MU spktimes
if ~isempty(v.spktimes)
    mutimes = cellfun(@(x) x / fsSpk, v.spktimes, 'uni', false);
else
    mutimes = [];
end
spkTimes = v.spikes.times;

% Run Analysis
rippSpks = ripp_spks(rippTimes, spkTimes, ...
    'basepath', basepath, ...
    'muTimes', mutimes, ...
    'flgPlot', flgPlot, ...
    'flgSave', flgSave);


%% ========================================================================
%  SPK-LFP
%  ========================================================================
fprintf('Running spklfp coupling...\n');
fRange = [120 200];
sigFilt = filterLFP(sig, 'fs', fs, 'type', 'butter', 'dataOnly', true, ...
    'order', 3, 'passband', fRange, 'graphics', false);

ripp.spkLfp = spklfp_calc('basepath', basepath, 'lfpTimes', ripp.times, ...
    'ch', ripp.info.rippCh, 'fRange', fRange, 'sig', sigFilt, ...
    'flgSave', false, 'flgPlot', false, 'flgStl', false);

if flgPlot
    spklfp_plot(ripp.spkLfp);
end

% Final Save (Core struct updates like spkLfp)
if flgSave
    save(fullfile(basepath, [basename, '.ripp.mat']), 'ripp', '-v7.3');
end

fprintf('--- Ripple Analysis Pipeline Completed for %s ---\n', basename);

end
