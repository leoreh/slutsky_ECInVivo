function ripp = ripp_wrapper(varargin)
% RIPP_WRAPPER Master function to run the complete ripple analysis pipeline.
%
% SUMMARY:
% This function serves as a wrapper to execute the core ripple analysis
% sequence: detection, state analysis, and spike analysis, along with their
% corresponding plotting functions. It allows for flexible parameter setting
% for each step through varargin.
%
% Pipeline Steps:
%   1. ripp_detect:  Detects ripples based on LFP, optional EMG, and thresholds.
%                    Can optionally refine detection using initial results.
%   2. ripp_plot:    (Optional) Generates summary plots for ripple detection.
%   3. ripp_states:  Analyzes ripple occurrence within vigilance states.
%                    Generates plots of state-specific rates.
%   4. ripp_spks:    Analyzes multi-unit (MU) and single-unit (SU) spiking
%                    relative to ripple events and matched control periods.
%                    Can limit statistical analysis to a specific state.
%   5. ripp_plotSpks:(Optional) Generates summary plots for spike analysis.
%
% INPUT (Optional Key-Value Pairs):
%   basepath        Path to recording session directory {pwd}.
%   rippCh          Channel(s) for ripple detection {session default}.
%   thr             Initial detection thresholds (see ripp_detect)
%                   {[1.5, 2.5, 2, 200, 100]}.
%   limDur          Initial duration limits (see ripp_detect)
%                   {[20, 300, 20, 10]}.
%   passband        Filter passband [Hz] {[100 300]}.
%   detectAlt       Detection method (see ripp_detect) {3}.
%   limState        Vigilance state index to limit ripp_spks statistics {[4]}.
%   flgRefine       Logical flag to run ripp_detect twice, using suggested
%                   parameters from the first run {true}.
%   flgPlot     Master flag to enable/disable ALL plotting {true}.
%   flgSave      Master flag to save/update ripp.mat at each step {true}.
%   flgSaveFig      Master flag to save ALL generated figures {true}.
%   bit2uv          Conversion factor from bits to microvolts {0.195}.
%
% OUTPUT:
%   ripp            Structure containing results from all analysis steps.
%                   Key top-level fields:
%                     .times, .peakTime, .dur, .peakPower, etc. (from detect)
%                     .info: Metadata about detection and analysis parameters.
%                     .maps: Peri-event LFP maps (from detect).
%                     .rate: Overall ripple rate (from detect).
%                     .acg: Autocorrelogram (from detect).
%                     .corr: Correlations between ripple features (from detect).
%                     .states: State-specific analysis (from states).
%                       .rate, .idx, .tstamps, etc.
%                     .spks: Spike analysis results (from spks).
%                       .mu: Multi-unit analysis (rippMap, ctrlMap).
%                       .su: Single-unit analysis (rippRates, ctrlRates,
%                            pVal, sigMod, spikeGain, rippMap, ctrlMap).
%                       .ctrlPeakTime, .ctrlTimes: Control event timings.
%
% DEPENDENCIES:
%   ripp_detect.m, ripp_plot.m, ripp_states.m, ripp_spks.m, ripp_plotSpks.m
%   All dependencies listed within those functions.
%
% HISTORY:
% Aug 2024 LH - Initial wrapper function creation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'rippCh', [], @isnumeric);
addParameter(p, 'thr', [1, 2.5, 1.2, 200, 100], @(x) isnumeric(x) && numel(x)==5);
addParameter(p, 'limDur', [20, 300, 20, 10], @(x) isnumeric(x) && numel(x)==4);
addParameter(p, 'passband', [100 300], @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'detectAlt', 3, @(x) isnumeric(x) && isscalar(x) && ismember(x, [1, 2, 3]));
addParameter(p, 'limState', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'flgRefine', true, @islogical);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgSave', true, @islogical);
addParameter(p, 'flgSaveFig', true, @islogical);
addParameter(p, 'bit2uv', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));

parse(p, varargin{:});
basepath        = p.Results.basepath;
rippCh          = p.Results.rippCh;
thr             = p.Results.thr;
limDur          = p.Results.limDur;
passband        = p.Results.passband;
detectAlt       = p.Results.detectAlt;
limState        = p.Results.limState;
flgRefine       = p.Results.flgRefine;
flgPlot         = p.Results.flgPlot;
flgSave      = p.Results.flgSave;
flgSaveFig      = p.Results.flgSaveFig;
bit2uv          = p.Results.bit2uv;

% Set basepath
cd(basepath);
[~, basename] = fileparts(basepath);
fprintf('--- Starting Ripple Analysis Pipeline for %s ---\n', basename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data params from session file
sessionfile = fullfile(basepath, [basename, '.session.mat']);
load(sessionfile, 'session');
fs = session.extracellular.srLfp;
fsSpk = session.extracellular.sr;
nchans = session.extracellular.nChannels;

% Get ripple channels
if isempty(rippCh)
    if isfield(session.channelTags, 'Ripple')
        rippCh = session.channelTags.Ripple;
    else
        rippCh = 1;
    end
end

% Handle bit2uv logic 
if isempty(bit2uv)
    if round(fsSpk) == 24414
        bit2uv = 1;      % Default for TDT
    else
        bit2uv = 0.195;  % Default for Intan/OE
    end
end


% EMG
load([basename, '.sleep_sig.mat'], 'emg');

% LFP
fname = fullfile(basepath, [basename, '.lfp']);
sig = binary_load(fname, 'duration', Inf,...
    'fs', fs, 'nCh', nchans, 'start', 0,...
    'ch', rippCh, 'downsample', 1, 'bit2uv', bit2uv);

if size(sig, 2) > 1
    sig = mean(sig, 2);
end

% Validate EMG length 
if length(emg) ~= length(sig)
     error('Provided EMG length (%d) does not match LFP length (%d)', length(emg), length(sig));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIPPLE DETECTION (+ OPTIONAL REFINEMENT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Running ripp_detect...\n');

if flgRefine
    winIdx = 1 : 3 * 60 * 60 * fs;
    thrTmp = [1, 2.5, 1.2, 200, 50];
    limDurTmp = [20, 300, 20, 10];
    refineFactor = 0.8;
    passband = [100 300];

    ripp = ripp_detect('basepath', basepath, 'sig', sig(winIdx), 'emg', emg(winIdx),...
        'rippCh', [], 'recWin', [0 Inf], 'thr', thrTmp, 'passband', passband,...
        'limDur', limDurTmp, 'detectAlt', detectAlt, 'bit2uv', bit2uv, ...
        'flgPlot', false, 'flgSave', false);

    % Use suggested parameters, adjusting thresholds
    thr = round(ripp.info.thrData * refineFactor, 1);
    thr(5) = 50; 
    limDur = ripp.info.limDurData;
end

ripp = ripp_detect('basepath', basepath, 'sig', sig, 'emg', emg,...
    'rippCh', [], 'recWin', [0 Inf], 'thr', thr, 'passband', passband,...
    'limDur', limDur, 'detectAlt', detectAlt, 'bit2uv', bit2uv, ...
    'flgPlot', false, 'flgSave', flgSave);

% update channel info
ripp.info.rippCh = rippCh;

% Plot detection results
if flgPlot
    fprintf('Running ripp_plot...\n');
    ripp_plot(ripp, 'basepath', basepath, 'flgSaveFig', flgSaveFig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIPPLE STATE ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Running ripp_states...\n');
ripp = ripp_states(ripp, 'basepath', basepath, 'flgPlot', flgPlot, ...
    'flgSave', flgSave, 'flgSaveFig', flgSaveFig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIPPLE SPIKE ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Running ripp_spks...\n');
ripp = ripp_spks(ripp, 'basepath', basepath, 'limState', limState,...
    'flgPlot', false, 'flgSave', flgSave);

% Plot spike analysis
if flgPlot
    ripp_plotSpks(ripp, 'basepath', basepath, 'flgSaveFig', flgSaveFig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPIKE-LFP COUPLING ANALYSIS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Running spklfp coupling...\n');

% filter lfp
fRange = [120 200];
sig = filterLFP(sig, 'fs', fs, 'type', 'butter', 'dataOnly', true,...
    'order', 3, 'passband', fRange, 'graphics', false);

ripp.spkLfp = spklfp_calc('basepath', basepath, 'lfpTimes', ripp.times,...
    'ch', ripp.info.rippCh, 'fRange', fRange, 'sig', sig,... 
    'flgSave', false, 'flgPlot', false, 'flgStl', false);

if flgPlot
    spklfp_plot(ripp.spkLfp)
end

% Save the updated ripp structure 
if flgSave
    save(fullfile(basepath, [basename, '.ripp.mat']), 'ripp', '-v7.3');
end

fprintf('--- Ripple Analysis Pipeline Completed for %s ---\n', basename);

end % EOF ripp_wrapper

