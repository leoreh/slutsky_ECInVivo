function ss = accusleep_wrapper(varargin)

% wrapper for state classification via AccuSleep by Yang Dan
% paper: Barger et al., PlosOne, 2019
% git: https://github.com/zekebarger/AccuSleep
% documentation: doc AccuSleep_instructions
% allows for user defined ignore times obtained automatically from the eeg
% spectrogram or manually using the AccuSleep gui
%
% INPUT:
%   basepath    string. path to recording folder {pwd}
%   cleanRec    numeric. used to remove artifacts manually via the gui (1)
%               or automatically via the spectrogram {2}. if empty will
%               classify entire recording. 
%   SR          numeric. sampling frequency of EMG and EEG {512}
%   epochLen    numeric. length of epoch [s] {2.5}
%   minBoutLen  numeric. length of minimum episodes. if empty will be equal
%               to epochLen
%   recSystem   string. recording system, {'tdt'} or 'oe'
%   calfile     string. path to calibrationData file. if empty will search
%               in mouse folder.
%   lfpCh       numeric. channel number of eeg to load from lfp file. can
%               be a vector and then the channels will be averaged 
%   emgCh       numeric. channel number of eeg to load from lfp file. for
%               oe recording system
%   viaGui          logical. perform analysis manually via gui (true) or
%                   automatically via this script {false}
%   forceCalibrate  logical. create calibration matrix even if one already
%                   exists for this mouse {false}
%   inspectLabels   logical. manually review classification 
%   saveVar         logical. save ss var {true}
%   forceAnalyze    logical. reanalyze recordings even if ss struct
%                   exists (false)
%   forceLoad       logical. reload recordings even if mat exists
%
% DEPENDENCIES:
%   AccuSleep
%   IOSR.DSP.SINCFILTER     for low-pass filtering EEG data
%
% TO DO LIST:
%       # reduce time bin of rec cleaning
%       # uigetf for net (done)
%       # improve removing ignoreEpochs (slow and inaccurate)
%       # graphics
%       # implement forceLoad (done)
%       # filter before resampling to assure nyquist (done not tested)
%       # nstates input argument
%       # save all output to accuslip folder
%       # check if mouse calibration exists
%
% 06 feb 21 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'cleanRec', '[]', @isnumeric);
addOptional(p, 'SR', 512, @isnumeric);
addOptional(p, 'epochLen', 2.5, @isnumeric);
addOptional(p, 'minBoutLen', [], @isnumeric);
addOptional(p, 'recSystem', 'tdt', @ischar);
addOptional(p, 'calfile', []);
addOptional(p, 'lfpCh', [], @isnumeric);
addOptional(p, 'emgCh', [], @isnumeric);
addOptional(p, 'viaGui', [], @islogical);
addOptional(p, 'forceCalibrate', [], @islogical);
addOptional(p, 'inspectLabels', [], @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'forceAnalyze', false, @islogical);
addOptional(p, 'forceLoad', false, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
cleanRec        = p.Results.cleanRec;
SR              = p.Results.SR;
epochLen        = p.Results.epochLen;
minBoutLen      = p.Results.minBoutLen;
recSystem       = p.Results.recSystem;
calfile         = p.Results.calfile;
lfpCh           = p.Results.lfpCh;
emgCh           = p.Results.emgCh;
viaGui          = p.Results.viaGui;
forceCalibrate  = p.Results.forceCalibrate;
inspectLabels   = p.Results.inspectLabels;
saveVar         = p.Results.saveVar;
forceAnalyze    = p.Results.forceLoad;
forceLoad       = p.Results.forceAnalyze;

cd(basepath)
if isempty(minBoutLen)
    minBoutLen = epochLen;
end

% badEpochs = [];
% recSystem = 'tdt';
% SR = 512;            
% epochLen = 2.5;        
% lfpCh = 13 : 16;
% viaGui = false;
% cleanRec = 2;
% forceCalibrate = false;
% minBoutLen = epochLen;

% change to user input or from configuration file
nstates = 4;        
ss.labelNames{1} = 'REM';
ss.labelNames{2} = 'WAKE';
ss.labelNames{3} = 'NREM';
ss.labelNames{4} = 'DROWSY';

% session info
[~, basename] = fileparts(basepath);
load([basename, '.session.mat'])
basepath = session.general.basePath;
nchans = session.extracellular.nChannels;
recDur = session.general.duration;
fsLfp = session.extracellular.srLfp;

% files
mousepath = fileparts(basepath);
calfile = [basename, '.AccuSleep_calibration.mat'];
callabelsfile = [basename, '.AccuSleep_labelsCalibration.mat'];
labelsfile = [basename, '.AccuSleep_labels.mat'];
eegfile = [basename '.AccuSleep_EEG.mat'];
emgfile = [basename '.AccuSleep_EMG.mat'];
artifactsfile = [basename, '.AccuSleep_labelsArtifacts.mat'];
statesfile = [basename '.AccuSleep_states.mat'];

netfile = 'D:\Code\AccuSleep\trainedNetworks\trainedNetwork2,5secEpochs.mat';
netfile = 'D:\Code\slutskycode\extracellular in vivo\lfp\detectStates\AccuSleep\3states_2,5s_6hrLabels_net.mat';
netfile = 'D:\Code\slutskycode\extracellular in vivo\lfp\detectStates\AccuSleep\4states_2,5s_12hrLabels_net.mat';
if ~exist(netfile)
    [netfile, netpath] = uigetfile;
    netfile = [fullfile(netpath, netfile), '.mat'];
end
load(netfile, 'net')
    
% check if already analyzed
if exist(statesfile, 'file') && ~forceAnalyze
    load(statesfile, 'ss')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(emgfile, 'file') && exist(eegfile, 'file') && ~forceLoad
    load(eegfile, 'EEG')
    load(emgfile, 'EMG')
else
    switch recSystem
        case 'tdt'
            % emg
            emgInfo = dir([basename, '.EMG*.datInfo*']);
            load(emgInfo.name)
            fsEmg = datInfo.fs;
            emgname = [basename '.emg.dat'];
            emg_orig = double(bz_LoadBinary(emgname, 'duration', Inf,...
                'frequency', fsEmg, 'nchannels', 1, 'start', 0,...
                'channels', 1, 'downsample', 1));
            
            % lfp (same for oe system)
            lfpname = [basename, '.lfp'];
            eeg_orig = double(bz_LoadBinary(lfpname, 'duration', Inf,...
                'frequency', fsLfp, 'nchannels', nchans, 'start', 0,...
                'channels', lfpCh, 'downsample', 1));
            if size(eeg_orig, 2) > 1
                eeg_orig = mean(eeg_orig, 2);
            end
        case 'oe'
    end
    
    % find corresponsding timestamps for EMG and labels [s]
    tstamps_sig = [1 / SR : 1 / SR : recDur];
    
    % low-pass filter to assure nyquist.
    % note accusleep only uses spectrogram up to 50 Hz
    %     cf = 100;        % cutoff frequency
    %     fprintf('\nlow-pass filtering, cutoff = %d Hz\n', cf)
    %     import iosr.dsp.*
    %     filtRatio = cf / (fsLfp / 2);
    %     eeg_orig = iosr.dsp.sincFilter(eeg_orig, filtRatio);
    
    % remove DC component
    % fprintf('removing dc component\n')
    % EEG = rmDC(EEG, 'dim', 1);
    % EEG = rmDC(EEG, 'dim', 1);
    
    % resmaple
    EMG = [interp1([1 : length(emg_orig)] / fsEmg, emg_orig, tstamps_sig,...
        'pchip')]';
    EEG = [interp1([1 : length(eeg_orig)] / fsLfp, eeg_orig, tstamps_sig,...
        'pchip')]';
    
    % save files
    save(eegfile, 'EEG')
    save(emgfile, 'EMG')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean recording from artifacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% two alternatives - semi-automatically or manually (see below)
% -------------------------------------------------------------------------
% a different approach is to mark artifacts during the calibration process
% as undefined. than, must insert the calibration labels in place of the
% final results or use the gui and select the option only overwrite
% undefined epochs. however, this will not prevent the artifacts from
% influincing the spectrogram if it is normalized. 
% -------------------------------------------------------------------------
% currently uses a temporal resolution that corresponds to epoch length but
% should be improved to < 1 s. just need to find a way to restore deleted
% sections afterwards.

cleanRes = 2.5;       % size of bins for cleaning the data [s]
artifactThr = 2;      % threshold for removing sig in alt 2 [z scores]
if exist(artifactsfile) && ~forceLoad
    load(artifactsfile)
else
    if cleanRec == 1
        % ALT 1: manual mark bad times using AccuSleep gui (use NREM
        % as a bad epoch). when done save labels.
        AccuSleep_viewer(EEG, EMG, SR, cleanRes, [], artifactsfile)
        uiwait
        load(artifactsfile)
        ignoretimes(labels == 3) = 1;
        ignoretimes(labels == 4) = 0;
        
        % ALT 2: semi-automatically find bad epochs from spectrogram
    elseif cleanRec == 2
        [s, ~, ~] = createSpectrogram(standardizeSR(EEG, SR, 128), 128, cleanRes);
        ignoretimes = zeros(floor(recDur / cleanRes), 1);
        ignoretimes(zscore(mean(zscore(s), 2)) > artifactThr) = 1;
    else
        ignoretimes = [];
    end
end
save(artifactsfile, 'ignoretimes')

% interpolate badtimes to length of the signals
if ~isempty(ignoretimes)
rmtimes = [interp1([1 : length(ignoretimes)] * cleanRes, ignoretimes, tstamps_sig,...
    'linear', 'extrap')]';
else
    rmtimes = zeros(length(EEG), 1);
end
rmtimes = round(rmtimes);
EEG = EEG(~rmtimes);
EMG = EMG(~rmtimes);
tRemoved = sum(rmtimes) / SR; % total time removed from data [s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AccuSleep pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labels = ones(1, floor(length(EMG) / SR / epochLen)) * nstates + 1;
save(labelsfile, 'labels')

if viaGui
    % gui
    AccuSleep_GUI
    uiwait
else
    % get mouse calibration. if does not exist already, manualy score some
    % of the data and save the labels
    if ~exist(calfile, 'file') || forceCalibrate
        if exist(callabelsfile)
            load(callabelsfile)
        end
        AccuSleep_viewer(EEG, EMG, SR, epochLen, labels, callabelsfile);
        uiwait
        load(callabelsfile)
        calibrationData = createCalibrationData(standardizeSR(EEG, SR, 128),...
            standardizeSR(EMG, SR, 128), labels, 128, epochLen);
        save(calfile, 'calibrationData')
    else
        load(calfile)
    end
    save(callabelsfile, 'labels')
    labels_calibration = labels;
    
    % classify recording
    labels_orig = AccuSleep_classify(EEG, EMG, net, SR, epochLen, calibrationData, minBoutLen);
    labels = labels_orig;
    save(labelsfile, 'labels')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restore bad times to labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tstamps_labels = [0 : epochLen : recDur];
% labelsLen = floor(length(EMG_newFs) / SR / epochLen);
% x = [interp1([1 : length(badtimes)] / cleanRes, badtimes,...
%     [1 : labelsLen] * epochLen, 'linear', 'extrap')]';
% x = round(x);

ignoreEpochs = binary2epochs('vec', ignoretimes, 'minDur', [], 'maxDur', [],...
    'interDur', [], 'exclude', false); 
for i = 1 : size(ignoreEpochs, 1)
    labels = [labels(1 : ignoreEpochs(i, 1) - 1);...
        ones(diff(ignoreEpochs(i, :)), 1) * 4;...
    labels(ignoreEpochs(i, 1) : end)];

    labels_calibration = [labels_calibration(1 : ignoreEpochs(i, 1) - 1);...
        ones(diff(ignoreEpochs(i, :)), 1) * 4;...
    labels_calibration(ignoreEpochs(i, 1) : end)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% review classification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measure accuracy compared to calibration data. output is a mat nstates x
% nstates where each row is the gold standard state (in this case calibration
% labels) and each column is the 
stateEpochs = histcounts(labels_calibration);
for j = 1 : nstates
    statesAccuracy(j, :) = histcounts(labels(labels_calibration == j),...
        [1 : nstates + 1]) / stateEpochs(j) * 100;
    fprintf('accuracy for %s: %.1f\n', ss.labelNames{j}, statesAccuracy(j, j))    
end

% insert calibrated labels to final results. this is because
% AccuSleep_classify doesn't allow for 'only overwrite undefind' as in the
% gui
calidx = find(labels_calibration ~= nstates + 1);
labels(calidx) = labels_calibration(calidx);

if inspectLabels
    AccuSleep_viewer(EEG, EMG, SR, epochLen, labels, labelsfile)
    uiwait
    load(labelsfile)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create states structure and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert labels to state epochs. 
for i = 1 : nstates
    binaryVec = zeros(length(labels), 1);
    binaryVec(labels == i) = 1;
    stateEpisodes = binary2epochs('vec', binaryVec, 'minDur', [], 'maxDur', [],...
        'interDur', [], 'exclude', false); % these are given as indices are equivalent to seconds
    ss.stateEpochs{i} = stateEpisodes * epochLen;
end

ss.labels = labels;
ss.labels_orig = labels_orig;
ss.labels_calibration = labels_calibration;
ss.calibrationData = calibrationData;
ss.rmtimes = rmtimes;
ss.tRemoved = tRemoved;
ss.fs = SR;
ss.net = netfile; 
ss.accuracy = statesAccuracy;

if saveVar
    save(statesfile, 'ss')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


