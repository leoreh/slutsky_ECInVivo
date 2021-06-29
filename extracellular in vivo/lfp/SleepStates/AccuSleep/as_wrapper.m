function ss = as_wrapper(EEG, EMG, sigInfo, varargin)

% wrapper for state classification via AccuSleep by Yang Dan
% paper: Barger et al., PlosOne, 2019
% git: https://github.com/zekebarger/AccuSleep
% documentation: doc AccuSleep_instructions
%
% INPUT:
%   EMG             numeric. emg data (1 x n)
%   EEG             numeric. eeg data (1 x n)
%   info            struct. see as_prepSig.
%   basepath        string. path to recording folder {pwd}
%   fs              numeric. sampling frequency of EMG and EEG {512}
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
%   saveFig         logical. save figure {true}
%
% DEPENDENCIES:
%   AccuSleep (modified in slutskycode)
%   IOSR.DSP.SINCFILTER     for filtering data
%
% TO DO LIST:
%       # filter before resampling to assure nyquist (done)
%       # implement forceLoad (done)
%       # uigetf for net (done)
%       # implement cleanSig 
%       # graphics
%
% 06 feb 21 LH  updates:
% 19 apr 21     separated prepSig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'fs', [], @isnumeric);
addOptional(p, 'calfile', []);
addOptional(p, 'netfile', []);
addOptional(p, 'viaGui', [], @islogical);
addOptional(p, 'forceCalibrate', [], @islogical);
addOptional(p, 'inspectLabels', [], @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'forceAnalyze', false, @islogical);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveFig', true, @islogical);

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
saveFig         = p.Results.saveFig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constants
epochLen = 1;        
minBoutLen = epochLen;

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
    netfile = [fullfile(netpath, netfile), '.mat'];
end
load(netfile, 'net')    

% validate data
if length(EEG) ~= length(EMG)
    error('EEG and EMG must be the same length')
end

if isempty(fs)
    fs = sigInfo.fs;
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
            labels = ones(1, floor(length(EMG) / fs / epochLen)) * nstates + 1;
        end
        
        fprintf('\nlabel some data for calibration, then press save.\n')
        AccuSleep_viewer(EEG, EMG, fs, epochLen, labels, callabelsfile);
        uiwait
        load(callabelsfile)
        
        % ignore bin state     
        labels_calibration = labels;
        labels(labels > nstates - 1) = nstates;          
        
        % calibration matrix
        fprintf('creating calbiration matrix... ')
        calibrationData = createCalibrationData(standardizeSR(EEG, fs, 128),...
            standardizeSR(EMG, fs, 128), labels, 128, epochLen);
        save(calfile, 'calibrationData')
        fprintf('done.\n')
    else
        load(calfile)
        load(callabelsfile)
        labels_calibration = labels;
    end
    
    % classify recording
    fprintf('classifying... ')
    [labels_net, netScores] = AccuSleep_classify(standardizeSR(EEG, fs, 128),...
        standardizeSR(EMG, fs, 128), net, 128, epochLen, calibrationData, minBoutLen);
    labels = labels_net;
    save(labelsfile, 'labels')
    fprintf('done.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    % confusion matrix (relative to clabiration data)
    labels1 = labels_calibration;
    labels2 = labels_net;
    [ss.netPrecision, ss.netRecall] = as_cm(labels1, labels2);

    % stateSeparation
    as_stateSeparation(EEG, EMG, labels)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finalize and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% insert calibrated labels to final results. this is because
% AccuSleep_classify doesn't allow for 'only overwrite undefind' as in the
% gui
calidx = find(labels_calibration ~= nstates + 1);
labels(calidx) = labels_calibration(calidx);

if inspectLabels
    AccuSleep_viewer(EEG, EMG, fs, epochLen, labels, labelsfile)
    uiwait
    load(labelsfile)
end

% convert labels to state epochs. 
for istate = 1 : nstates
    binaryVec = zeros(length(labels), 1);
    binaryVec(labels == istate) = 1;
    stateEpochs = binary2epochs('vec', binaryVec, 'minDur', [], 'maxDur', [],...
        'interDur', [], 'exclude', false); % these are given as indices and are equivalent to seconds
    ss.stateEpochs{istate} = stateEpochs * epochLen;
    ss.epLen{istate} = [diff(ss.stateEpochs{istate}')]';
end

ss.info = sigInfo;
ss.info.net = netfile;
ss.info.calibrationData = calibrationData;
ss.info.analysisDate = datetime;
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

