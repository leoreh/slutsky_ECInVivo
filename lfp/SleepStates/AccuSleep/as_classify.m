function ss = as_classify(sSig, varargin)

% state classification via AccuSleep by Yang Dan paper: Barger et al.,
% PlosOne, 2019 git: https://github.com/zekebarger/AccuSleep documentation:
% doc AccuSleep_instructions. as_prepSig must be called beforehand.
% function is basically a wrapper for AccuSleep_classify. 
%
% INPUT:
%   sSig            struct. see as_prepSig.m
%   basepath        string. path to recording folder {pwd}
%   netfile         string. path to network file. 
%   calData         numeric. calibration data. 
%   inspectLabels   logical. manually review classification 
%   saveVar         logical. save ss var {true}
%   forceA          logical. reanalyze recordings even if ss struct
%                   exists (false)
%   graphics        logical. plot confusion chart and state separation {true}
%
% DEPENDENCIES:
%   AccuSleep (modified in slutskycode)
%   IOSR.DSP.SINCFILTER     for filtering data
%   binary2epochs
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
% 06 jan 22     cleaned and implemented sSig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'netfile', []);
addOptional(p, 'calData', [], @isnumeric);
addOptional(p, 'inspectLabels', [], @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'forceA', false, @islogical);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
netfile         = p.Results.netfile;
calData         = p.Results.calData;
inspectLabels   = p.Results.inspectLabels;
saveVar         = p.Results.saveVar;
forceA          = p.Results.forceA;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load config data
cfg = as_loadConfig();
ss.info.names = cfg.names;
ss.info.colors = cfg.colors;
ss.info.epochLen = cfg.epochLen;
ss.info.minBoutLen = cfg.minBoutLen;
nstates = cfg.nstates;

% arrange files 
cd(basepath)
[~, basename] = fileparts(basepath);
manlabelsfile = [basename, '.sleep_labelsMan.mat'];
labelsfile = [basename, '.sleep_labels.mat'];
statesfile = [basename '.sleep_states.mat'];

% check if already analyzed
if exist(statesfile, 'file') && ~forceA
    fprintf('%s already exist. loading...', statesfile)
    load(statesfile, 'ss')
    return
end

% network file
if isempty(netfile)
    netpath = 'D:\Code\slutsky_ECInVivo\lfp\SleepStates\AccuSleep\trainedNetworks\';
    netfiles = dir([netpath, '/net*']);
    netfiles = sort({netfiles.name});
    netfile = netfiles{end};
    netfile = [fullfile(netpath, netfile)];
end
if ~exist(netfile, 'file')
    [netfile, netpath] = uigetfile('', 'Please select network file');
    netfile = [fullfile(netpath, netfile)];
end
load(netfile, 'net', 'netInfo')    

% validate data
if length(sSig.eeg) ~= length(sSig.emg)
    error('EEG and EMG must be the same length')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AccuSleep pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load manual labels if exists
if exist(manlabelsfile, 'file')
    load(manlabelsfile, 'labels')

    % ignore bin state
    labels_man = labels;
    labels(labels > nstates) = nstates + 1;

else
    labels_man = ones(1, length(sSig.spec_tstamps)) * nstates + 2;
end

% get mouse calibration or manualy score some of the data and save the labels
if isempty(calData)
    if ~exist(manlabelsfile, 'file')
        fprintf('\nlabel some data for calibration, then press save.\n')
        AccuSleep_viewer(sSig, [], manlabelsfile)
        uiwait
        load(manlabelsfile, 'labels')
        labels_man = labels;
        labels(labels > nstates) = nstates + 1;
    end
  
    % calibration matrix
    fprintf('\ncreating calbiration matrix... ')
    calData =...
        createCalibrationData(sSig.spec, sSig.spec_freq, sSig.emg_rms, labels);
    fprintf('done.\n')
end

% classify recording
fprintf('classifying... ')
[labels_net, netScores] = AccuSleep_classify(sSig.spec,...
    sSig.spec_freq, sSig.emg_rms, calData, net);

% insert munual labels to final results. this is because
% AccuSleep_classify doesn't allow for 'only overwrite undefind' as in the
% gui
labels = labels_net;
manIdx = find(labels_man < nstates + 1);
labels(manIdx) = labels_man(manIdx);
fprintf('done.\n')

% remove labels of uncertainty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert labels to state epochs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minDur = [10, 5, 5, 10, 5, 5];
interDur = 4;

[ss.stateEpochs, epochStats] = as_epochs('labels', labels,...
    'minDur', minDur, 'interDur', interDur);
ss.epLen = epochStats.epLen;
ss.nepochs = epochStats.nepochs;
ss.totDur = epochStats.totDur;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finalize and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ss.info.sSig = sSig.info;
ss.info.net = netfile;
ss.info.calibrationData = calData;
ss.info.minDur = minDur;
ss.info.interDur = interDur;
ss.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
ss.labels = labels;
ss.labels_net = labels_net;
ss.labels_man = labels_man;
ss.netScores = netScores;

if saveVar
    save(statesfile, 'ss')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    % confusion matrix (relative to manual labels)
    if any(labels_man ~= nstates + 2)
        [ss.netPrecision, ss.netRecall] = as_cm(labels_man, labels_net, netScores);
    end
    
    % stateSeparation
    as_stateSeparation(sSig, ss)
    
    % duration in timebins
    sstates = find(matches(cfg.names, {'WAKE', 'NREM', 'REM'}));
    as_plotZT('nbins', 4, 'sstates', sstates, 'ss', ss)
    
end

if inspectLabels
    AccuSleep_viewer(sSig, labels, labelsfile)
    uiwait
    if exist(labelsfile, 'file')
        load(labelsfile, 'labels')
    end
end

return

% EOF

% -------------------------------------------------------------------------
% subroutine for combining lables from two different classifications

nlabels = length(ss.labels);

idx1 = [1 : 12 * 60 * 60];
idx2 = [12 * 60 * 60 + 1 : nlabels];

labels = ones(nlabels, 1);
labels(idx1) = ss.labels(idx1);
labels(idx2) = ss.labels(idx2);

