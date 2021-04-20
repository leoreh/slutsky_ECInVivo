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
configfile = 'D:\Code\AccuSleep\AS_config.mat';
if ~exist(configfile, 'file')
    scriptfile = mfilename('fullpath');
    scriptpath = fileparts(scriptfile);
    configfile(fullfile(scriptpath, 'AS_config.mat'));
    
    if ~exist(configfile, 'file')
        [configfile, configpath] = uigetfile('', 'Please select  the configuration file');
        configfile = [fullfile(configpath, configfile), '.mat'];
    end
end
load(configfile)
ss.labelNames = cfg_names;
nstates = length(ss.labelNames);

% network file
if isempty(netfile)
    netfile = 'D:\Code\AccuSleep\trainedNetworks\trainedNetwork2,5secEpochs.mat';
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
labels = ones(1, floor(length(EMG) / fs / epochLen)) * nstates + 1;

if viaGui
    AccuSleep_GUI       % open gui
    uiwait
else
    % get mouse calibration. if does not exist already, manualy score some
    % of the data and save the labels
    if ~exist(calfile, 'file') || forceCalibrate
        if exist(callabelsfile)
            load(callabelsfile)
        end
        AccuSleep_viewer(EEG, EMG, fs, epochLen, [], callabelsfile);
        uiwait
        load(callabelsfile)
        calibrationData = createCalibrationData(standardizeSR(EEG, fs, 128),...
            standardizeSR(EMG, fs, 128), labels, 128, epochLen);
        save(calfile, 'calibrationData')
    else
        load(calfile)
    end
    save(callabelsfile, 'labels')
    labels_calibration = labels;
    
    % classify recording
    [labels_net, netScores] = AccuSleep_classify(EEG, EMG, net, fs, epochLen, calibrationData, minBoutLen);
    labels = labels_net;
    save(labelsfile, 'labels')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% confusion matrix (relative to clabiration data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cm = confusionmat(labels_calibration, labels);

% function handles. precision: when predicts yes, how often is it correct?
% TP / (TP + FP). recall: when actually yes, how often does it predict yes?
% TP / (TP + FN)
precision = @(cm) diag(cm)./sum(cm, 2);
recall = @(confusionMat) diag(confusionMat) ./ sum(confusionMat, 1)';

% calc
ss.netPrecision = precision(cm);
ss.netRecall = recall(cm);

% plot confusion chart
if graphics
    % ALT 1:
    % plotconfusion(categorical(labels_calibration), categorical(labels_net))
    
    % ALT 2:
    fh = figure;
    fh.Position = [500 200 900 700];
    confusionchart(cm, cfg_names, 'ColumnSummary',...
        'column-normalized', 'RowSummary', 'row-normalized',...
        'title', 'State Classification Confusion Matrix', 'Normalization',...
        'total-normalized')
    
    if saveFig
        figpath = fullfile('graphics', 'sleepState');
        mkdir(figpath)
        figname = fullfile(figpath, sprintf('%s_stateConfusionMat', basename));
        export_fig(figname)
    end
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
for i = 1 : nstates
    binaryVec = zeros(length(labels), 1);
    binaryVec(labels == i) = 1;
    stateEpisodes = binary2epochs('vec', binaryVec, 'minDur', [], 'maxDur', [],...
        'interDur', [], 'exclude', false); % these are given as indices and are equivalent to seconds
    ss.stateEpochs{i} = stateEpisodes * epochLen;
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

end

% EOF

% -------------------------------------------------------------------------
% vars for direct running
calfile = [];
viaGui = false;
forceCalibrate = true;
inspectLabels = true;
saveVar = true;
forceAnalyze = true;
fs = 1000;



