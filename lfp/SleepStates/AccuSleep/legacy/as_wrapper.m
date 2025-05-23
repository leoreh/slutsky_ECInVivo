function ss = as_wrapper(EEG, EMG, sigInfo, varargin)

% wrapper for state classification via AccuSleep by Yang Dan
% paper: Barger et al., PlosOne, 2019
% git: https://github.com/zekebarger/AccuSleep
% documentation: doc AccuSleep_instructions
% ignores last state (assumes bin)
%
% INPUT:
%   EMG             numeric. emg data (1 x n)
%   EEG             numeric. eeg data (1 x n)
%   info            struct. see as_prepSig.
%   basepath        string. path to recording folder {pwd}
%   fs              numeric. sampling frequency of EMG and EEG 
%   calfile         string. path to calibrationData file. if empty will search
%                   in mouse folder.
%   netfile         string. path to network file. 
%   viaGui          logical. perform analysis manually via gui (true) or
%                   automatically via this script {false}
%   forceCalibrate  logical. create calibration matrix even if one already
%                   exists for this mouse {false}
%   inspectLabels   logical. manually review classification 
%   saveVar         logical. save ss var {true}
%   forceAnalyze    logical. reanalyze recordings even if ss struct
%                   exists (false)
%   graphics        logical. plot confusion chart and state separation {true}
%
% DEPENDENCIES:
%   AccuSleep (modified in slutskycode)
%   IOSR.DSP.SINCFILTER     for filtering data
%   binary2bouts
%
% TO DO LIST:
%       # filter before resampling to assure nyquist (done)
%       # implement forceLoad (done)
%       # uigetf for net (done)
%       # graphics (done)
%       # implement tsa filter
%       # batch processing
%       # load sig instead of input (memory)
%
% 06 feb 21 LH  updates:
% 19 apr 21     separated prepSig
% 29 nov 21     add minDur and interDur

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'fs', 1250, @isnumeric);
addOptional(p, 'calfile', []);
addOptional(p, 'netfile', []);
addOptional(p, 'viaGui', [], @islogical);
addOptional(p, 'forceCalibrate', [], @islogical);
addOptional(p, 'inspectLabels', [], @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'forceAnalyze', false, @islogical);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
fs              = p.Results.fs;
netfile         = p.Results.netfile;
calfile         = p.Results.calfile;
viaGui          = p.Results.viaGui;
forceCalibrate  = p.Results.forceCalibrate;
inspectLabels   = p.Results.inspectLabels;
saveVar         = p.Results.saveVar;
forceAnalyze    = p.Results.forceAnalyze;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constants
minEpochLen = 1;        
minBoutLen = minEpochLen;

% arrange files 
cd(basepath)
mousepath = fileparts(basepath);
[~, basename] = fileparts(basepath);
callabelsfile = [basename, '.AccuSleep_labelsCalibration.mat'];
manlabelsfile = [basename, '.AccuSleep_labelsMan.mat'];
labelsfile = [basename, '.AccuSleep_labels.mat'];
statesfile = [basename '.AccuSleep_states.mat'];
if isempty(calfile)
    calfile = [basename, '.AccuSleep_calibration.mat'];
end

% check if already analyzed
if exist(statesfile, 'file') && ~forceAnalyze
    fprintf('%s already exist. loading...', statesfile)
    load(statesfile, 'ss')
    return
end

% get params from configuration file
[cfg_colors, cfg_names, cfg_weights] = as_loadConfig([]);
ss.labelNames = cfg_names;
nstates = length(ss.labelNames);

% network file
if isempty(netfile)
    netpath = 'D:\Code\slutskycode\extracellular in vivo\lfp\SleepStates\AccuSleep\trainedNetworks';
    netfiles = dir([netpath, '/net*']);
    netfiles = sort({netfiles.name});
    netfile = netfiles{end};
end
if ~exist(netfile)
    [netfile, netpath] = uigetfile('', 'Please select network file');
    netfile = [fullfile(netpath, netfile)];
end
load(netfile, 'net')    

% validate data
if length(EEG) ~= length(EMG)
    error('EEG and EMG must be the same length')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AccuSleep pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if viaGui
    AccuSleep_GUI       % open gui
    uiwait
else
    % get mouse calibration. if does not exist already, manualy score some
    % of the data and save the labels
    if ~exist(callabelsfile, 'file') || forceCalibrate
        if exist(manlabelsfile)
            load(manlabelsfile)
        elseif exist(callabelsfile)
            load(callabelsfile)
        else
            fprintf('\nlabel some data for calibration, then press save.\n')
            AccuSleep_viewer(EEG, EMG, fs, minEpochLen, [], callabelsfile);
            uiwait
            load(callabelsfile)
        end

        % ignore bin state     
        labels_calibration = labels;
        labels(labels > nstates - 1) = nstates;          
        
        % calibration matrix
        fprintf('creating calbiration matrix... ')
        calibrationData = createCalibrationData(EEG,...
            EMG, labels, fs, minEpochLen);
        save(calfile, 'calibrationData')
        fprintf('done.\n')
    else
        load(calfile)
        load(callabelsfile)
        labels_calibration = labels;
    end
    
    % classify recording
    fprintf('classifying... ')
    [labels_net, netScores] = AccuSleep_classify(EEG,...
        EMG, net, fs, minEpochLen, calibrationData, minBoutLen);
    
    % insert calibrated labels to final results. this is because
    % AccuSleep_classify doesn't allow for 'only overwrite undefind' as in the
    % gui
    labels = labels_net;
    calidx = find(labels_calibration ~= nstates + 1);
    labels(calidx) = labels_calibration(calidx);
    save(labelsfile, 'labels')
    fprintf('done.\n')
end

x = labels;
idx = 11 * 60 * 60 : length(labels);
x(idx) = labels(idx);
labels = x;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert labels to state bouts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minDur = [10, 5, 5, 10, 5, 5];  % threshold of minimum bout length
interDur = 4;                   % combine bouts separated by <= interDur

if length(minDur) == 1
    minDur = repamt(minDur, nstates - 1, 1);
elseif length(minDur) ~= nstates - 1
    error('minDur length is different than the number of states')
end
if length(interDur) == 1
    interDur = repmat(interDur, nstates - 1, 1);
elseif length(interDur) ~= nstates - 1
    error('interDur length is different than the number of states')
end

% create state bouts
for istate = 1 : nstates - 1
    binaryVec = zeros(length(labels), 1);
    binaryVec(labels == istate) = 1;
    boutTimes = binary2bouts('vec', binaryVec, 'minDur', minDur(istate), 'maxDur', [],...
        'interDur', interDur(istate), 'exclude', false);
    ss.boutTimes{istate} = boutTimes * minEpochLen; % convert indices to seconds
    ss.boutLen{istate} = [diff(ss.boutTimes{istate}')]';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    % confusion matrix (relative to calibration data)
    labels1 = labels_calibration;
    labels2 = labels_net;
    [ss.netPrecision, ss.netRecall] = as_cm(labels1, labels2, netScores);

    % stateSeparation
    as_stateSeparation(EEG, EMG, labels, 'boutTimes', ss.boutTimes,...
        'fs', fs)
end

if inspectLabels
    AccuSleep_viewer(EEG, EMG, fs, minEpochLen, labels, labelsfile)
    uiwait
    load(labelsfile)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finalize and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ss.info = sigInfo;
ss.info.net = netfile;
ss.info.calibrationData = calibrationData;
ss.info.analysisDate = datetime;
ss.info.minEpochLen = minEpochLen;
ss.info.minBoutLen = minBoutLen;
ss.info.minDur = minDur;
ss.info.interDur = interDur;
ss.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
ss.labels = labels;
ss.labels_net = labels_net;
ss.labels_calibration = labels_calibration;
ss.netScores = netScores;

if saveVar
    save(statesfile, 'ss')
end

return

% EOF

% -------------------------------------------------------------------------
% vars for direct running
calfile = [];
viaGui = false;
forceCalibrate = true;
inspectLabels = true;
saveVar = true;
forceAnalyze = true;
fs = 1250;
netfile = [];
sigInfo = [];

