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
addParameter(p, 'thr', [1, 2.5, 1.2, 200, 50], @isnumeric);
addParameter(p, 'limDur', [15, 300, 30, 10], @isnumeric);
addParameter(p, 'win', [0, Inf], @isnumeric);
addParameter(p, 'passband', [80 250], @isnumeric);
addParameter(p, 'detectMet', 4, @isnumeric);
addParameter(p, 'limState', [], @isnumeric);
addParameter(p, 'flgGui', false, @islogical);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgSave', false, @islogical);
addParameter(p, 'flgSaveFig', false, @islogical);
addParameter(p, 'bit2uv', [], @isnumeric);
addParameter(p, 'zMet', 'nrem', @ischar);
parse(p, varargin{:});

basepath = p.Results.basepath;
rippCh = p.Results.rippCh;
thr = p.Results.thr;
limDur = p.Results.limDur;
passband = p.Results.passband;
win = p.Results.win;
detectMet = p.Results.detectMet;
flgPlot = p.Results.flgPlot;
flgSave = p.Results.flgSave;
flgSaveFig = p.Results.flgSaveFig;
bit2uv = p.Results.bit2uv;
limState = p.Results.limState;
flgGui = p.Results.flgGui;
zMet = p.Results.zMet;

%% ========================================================================
%  SETUP
%  ========================================================================

cd(basepath);
[~, basename] = fileparts(basepath);
fprintf('--- Ripple Analysis for %s ---\n', basename);

% Load vars
vars = {'session', 'spikes', 'spktimes', 'sleep_states'};
v = basepaths2vars('basepaths', {basepath}, 'vars', vars);

% Session params
fs = v.session.extracellular.srLfp;
fsSpk = v.session.extracellular.sr;
nchans = v.session.extracellular.nChannels;

if isinf(win(2))
    sigDur = Inf;
else
    sigDur = win(2) - win(1);
end

% Prep MU spktimes
if ~isempty(v.spktimes)
    muTimes = cellfun(@(x) x / fsSpk, v.spktimes, 'uni', false);
    muTimes = cellfun(@(x) x - win(1), muTimes, 'Uni', false);
    muTimes = cellfun(@(x) x(x >= 0 & x <= sigDur), muTimes, 'Uni', false);
    muTimes = {sort(vertcat(muTimes{:}))};
else
    muTimes = [];
end

% Prep SU spktimes
spkTimes = v.spikes.times;
spkTimes = cellfun(@(x) x - win(1), spkTimes, 'Uni', false);
spkTimes = cellfun(@(x) x(x >= 0 & x <= sigDur), spkTimes, 'Uni', false);

% Prep boutTimes (for state analysis)
try
    boutTimes = v.ss.bouts.times;
    boutTimes = cellfun(@(x) x - win(1), boutTimes, 'UniformOutput', false);
    boutTimes = cellfun(@(x) x(x(:,2)>0 & x(:,1)<sigDur, :), boutTimes, 'UniformOutput', false);
    nremTimes = boutTimes{4};
catch
    boutTimes = [];
    nremTimes = [];
end

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

%% ========================================================================
%  LOAD SIGNALS
%  ========================================================================

% LFP
fname = fullfile(basepath, [basename, '.lfp']);

lfp = binary_load(fname, 'duration', sigDur, 'fs', fs, 'nCh', nchans,...
    'start', win(1), 'ch', rippCh, 'downsample', 1, 'bit2uv', bit2uv);

if size(lfp, 2) > 1
    lfp = mean(lfp, 2);
end

% EMG
load([basename, '.sleep_sig.mat'], 'emg');

s1 = round(win(1) * fs) + 1;
if isinf(win(2))
    emg = emg(s1:end);
else
    s2 = min(length(emg), round(win(2) * fs));
    emg = emg(s1:s2);
end

if length(emg) ~= length(lfp)
    error('EMG and LFP do not fit')
end

% Prepare Signals
rippSig = ripp_sigPrep(lfp, fs, ...
    'detectMet', detectMet, 'passband', passband, ...
    'zMet', zMet, 'nremTimes', nremTimes);


%% ========================================================================
%  DETECT
%  ========================================================================

% Detect
ripp = ripp_detect(rippSig, fs, ...
    'emg', emg, ...
    'limDur', limDur, ...
    'basepath', basepath, ...
    'thr', thr, ...
    'flgPlot', false, ...
    'flgSave', flgSave);

% Get non-ripple intervals of matched duration.
recDur = length(lfp) / fs;
ctrlTimes = ripp_ctrlTimes(ripp.times, recDur);


if flgGui
    ripp_gui(ripp.times, ripp.peakTime, rippSig, muTimes, fs, thr, emg, ...
        'basepath', basepath);
end



%% ========================================================================
%  STATES
%  ========================================================================

rippStates = ripp_states(ripp.times, ripp.peakTime, boutTimes, ...
    'basepath', basepath, ...
    'flgPlot', flgPlot, ...
    'flgSave', flgSave, ...
    'flgSaveFig', flgSaveFig);

% Filter Ripples by State if requested
rippTimes = ripp.times;
peakTime = ripp.peakTime;
% if ~isempty(limState)
%     stateIdx = ripp.states.idx(:, limState);
%     rippTimes = ripp.times(stateIdx, :);
%     peakTime = ripp.peakTime(stateIdx);
%     fprintf('Limiting analysis to State %d\n', limState);
% end


%% ========================================================================
%  SPIKES
%  ========================================================================

% Run Analysis
rippSpks = ripp_spks(rippTimes, spkTimes, peakTime, ctrlTimes, ...
    'muTimes', muTimes, ...
    'basepath', basepath, ...
    'flgPlot', flgPlot, ...
    'flgSave', flgSave);



%% ========================================================================
%  SPK-LFP
%  ========================================================================
fprintf('Running spklfp coupling...\n');
fRange = [80 250];
% sig is not available anymore as vector, extract from rippSig or re-filter?
% rippSig.filt is [100 300]. Here we want [120 200].
% So we re-filter raw lfp.
sig = filterLFP(lfp, 'fs', fs, 'type', 'butter', 'dataOnly', true, ...
    'order', 3, 'passband', fRange, 'graphics', false);

spkLfp = spklfp_calc('basepath', basepath, ...
    'lfpTimes', ripp.times, ...
    'sig', rippSig.filt, ...
    'flgSave', false, ...
    'flgPlot', flgPlot, ...
    'flgStl', false);

spkLfpFile = fullfile(basepath, [basename, '.rippSpkLfp.mat']);
save(spkLfpFile, 'spkLfp', '-v7.3');


%% ========================================================================
%  Finalize and Save
%  ========================================================================

% NeuroScope Files
% Revert times to absolute for saving/Neuroscope
ripp.times = ripp.times + win(1);
ripp.peakTime = ripp.peakTime + win(1);

% NeuroScope Files
if flgSave

    % Convert to samples (NeuroScope uses acquisition rate fsSpk)
    rippSamps = round(ripp.times * fsSpk);
    peakSamps = round(ripp.peakTime * fsSpk);
    ripp2ns(rippSamps, peakSamps, 'basepath', basepath);
end

fprintf('--- Ripple Analysis Pipeline Completed for %s ---\n', basename);

% Put together all ripp structs
ripp.states = rippStates;
ripp.spks = rippSpks;
ripp.spkLfp = spkLfp;


end
