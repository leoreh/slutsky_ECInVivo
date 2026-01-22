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
addParameter(p, 'flgGui', false, @islogical);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgSave', false, @islogical);
addParameter(p, 'flgQA', true, @islogical);
addParameter(p, 'flgSaveFig', false, @islogical);
addParameter(p, 'bit2uv', [], @isnumeric);
addParameter(p, 'zMet', 'nrem', @ischar);
addParameter(p, 'mapDur', [-0.075 0.075], @isnumeric);
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
flgQA = p.Results.flgQA;
flgSaveFig = p.Results.flgSaveFig;
bit2uv = p.Results.bit2uv;
flgGui = p.Results.flgGui;
zMet = p.Results.zMet;
mapDur = p.Results.mapDur;

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
muTimes = cellfun(@(x) x / fsSpk, v.spktimes, 'uni', false);
muTimes = cellfun(@(x) x - win(1), muTimes, 'Uni', false);
muTimes = cellfun(@(x) x(x >= 0 & x <= sigDur), muTimes, 'Uni', false);
muTimes = {sort(vertcat(muTimes{:}))};

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
%  PRELIMINARY DETECTION
%  ========================================================================

% Detect
ripp = ripp_times(rippSig, fs, ...
    'thr', thr, ...
    'limDur', limDur);

% Get non-ripple intervals of matched duration.
ripp.ctrlTimes = ripp_ctrlTimes(ripp.times);


%% ========================================================================
%  PARAMS & MAPS
%  ========================================================================

ripp = ripp_params(rippSig, ripp);
rippMaps = ripp_maps(rippSig, ripp.peakTime, fs, 'mapDur', mapDur);
 
%% ========================================================================
%  STATES
%  ========================================================================

ripp.state = ripp_states(ripp.times, ripp.peakTime, boutTimes, ...
    'basepath', basepath, ...
    'flgPlot', false, ...
    'flgSave', false);

if flgSave
    rippFile = fullfile(basepath, [basename, '.ripp.mat']);
    save(rippFile, 'ripp', '-v7.3');
end


%% ========================================================================
%  SPIKE MAPS
%  ========================================================================

fprintf('Creating Spike PETH...\n');

rippPeth.su = ripp_spkPeth(spkTimes, ripp.peakTime, ripp.ctrlTimes, ...
    'flgSave', false, ...
    'mapDur', mapDur);
rippPeth.mu = ripp_spkPeth(muTimes, ripp.peakTime, ripp.ctrlTimes, ...
    'flgSave', false, ...
    'mapDur', mapDur);

if flgSave
    spkPethFile = fullfile(basepath, [basename, '.rippPeth.mat']);
    save(spkPethFile, 'rippPeth', '-v7.3');
end


%% ========================================================================
%  QUALITY ASSURANCE
%  ========================================================================

% Calculate Instantaneous Rates (Events x Units)
% 'binsize', Inf -> Returns one rate per event
rippRates = times2rate(muTimes, 'winCalc', ripp.times, 'binsize', Inf);
ctrlRates = times2rate(muTimes, 'winCalc', ripp.ctrlTimes, 'binsize', Inf);

% Calculate MUA spkGain (z-score relative to control) per ripple 
ripp.spkGain = [(rippRates - mean(ctrlRates, 'all', 'omitnan')) ./ ...
    std(ctrlRates, [], 'all', 'omitnan')]';

% Good ripples increase spiking activity
ripp.goodSpk = ripp.spkGain > 0;

if flgGui
    
    % MEAN TRACES
        
    % Build table
    tbl = struct2table(rmfield(ripp, {'info', 'times', 'ctrlTimes'}));
    
    % Add maps
    tblMap = struct2table(rmfield(rippMaps, {'tstamps'}));
    tblVars = tblMap.Properties.VariableNames;
    tblVars = strcat('t_', tblVars);
    tblMap.Properties.VariableNames = tblVars;
    tbl = [tbl, tblMap];
    
    % Add PETH
    tbl.suPETH = squeeze(mean(rippPeth.su.ripp, 1, 'omitnan'));
    tbl.muPETH = squeeze(rippPeth.mu.ripp);

    % Plot    
    tblGUI_xy(rippPeth.mu.tstamps, tbl, 'yVar', 't_z', 'grpVar', 'states');

    % PER TRACE
    ripp_gui(ripp.times, ripp.peakTime, rippSig, muTimes, fs, thr, emg, ...
        'basepath', basepath);
end

% Filter good ripples
if flgQA
    stateIdx = ripp.state == 'QWAKE' | ripp.state == 'LSLEEP' | ripp.state == 'NREM';
    goodIdx = stateIdx & ripp.goodSpk;
else
    goodIdx = true(length(ripp.state), 1);
end
ripp.goodIdx = goodIdx;

% Filter ripp (not saved, only for subsequent analysis)
fnames = fieldnames(ripp);
nRipp = length(goodIdx);
for iField = 1:length(fnames)
    fn = fnames{iField};
    
    % Check if first dimension matches nEvents
    if size(ripp.(fn), 1) == nRipp
        ripp.(fn) = ripp.(fn)(goodIdx, :);
    end
end
fprintf('Kept %d / %d events after QA.\n', sum(goodIdx), nRipp);

% Filter rippPeth (not saved, only for plotting)
fnames = fieldnames(rippPeth.su);
nRipp = length(goodIdx);
for iField = 1:length(fnames)
    fn = fnames{iField};
    
    % Check if first dimension matches nEvents
    if size(rippPeth.su.(fn), 2) == nRipp
        rippPeth.su.(fn) = rippPeth.su.(fn)(:, goodIdx, :);
    end

    % Check if first dimension matches nEvents
    if size(rippPeth.mu.(fn), 2) == nRipp
        rippPeth.mu.(fn) = rippPeth.mu.(fn)(:, goodIdx, :);
    end
end

% Recalculate states
ripp.state = ripp_states(ripp.times, ...
    ripp.peakTime, ...
    boutTimes, ...
    'basepath', basepath, ...
    'flgPlot', false, ...
    'flgSave', true);


%% ========================================================================
%  SPIKE STATS
%  ========================================================================

% SU
fprintf('Analyzing Spike Stats...\n');
rippSpks.su = ripp_spks(spkTimes, ...
    ripp.times, ...
    ripp.ctrlTimes, ...
    'basepath', basepath, ...
    'flgSave', false);

% MU
rippSpks.mu = ripp_spks(muTimes, ...
    ripp.times, ...
    ripp.ctrlTimes, ...
    'basepath', basepath, ...
    'flgSave', false);

if flgSave
    rippSpksFile = fullfile(basepath, [basename, '.rippSpks.mat']);
    save(rippSpksFile, 'rippSpks', '-v7.3');
end

% Plot spikes
if flgPlot 
    ripp_plotSpks(rippSpks, rippPeth, ...
        'basepath', basepath, ...
        'flgSaveFig', flgSaveFig);  
end


%% ========================================================================
%  SPK-LFP
%  ========================================================================
fprintf('Running spklfp coupling...\n');

spkLfp = spklfp_calc('basepath', basepath, ...
    'lfpTimes', ripp.times, ...
    'sig', rippSig.filt, ...
    'flgSave', false, ...
    'flgPlot', flgPlot, ...
    'flgStl', false);

spkLfpFile = fullfile(basepath, [basename, '.rippSpkLfp.mat']);
save(spkLfpFile, 'spkLfp', '-v7.3');


%% ========================================================================
%  NEUOROSCOPE
%  ========================================================================

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


%% ========================================================================
%  FINALIZE
%  ========================================================================

fprintf('--- Ripple Analysis Pipeline Completed for %s ---\n', basename);


end
