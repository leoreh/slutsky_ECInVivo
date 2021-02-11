% accusleep_simulation. trains a network on half a data set and test
% accuracy on the other half

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = 'D:\Data\lh81\lh81_210208_065300';
[~, basename] = fileparts(basepath);

labelsfile = [basename, '.AccuSleep_labels.mat'];
labelsmanfile = [basename, '.AccuSleep_labelsMan.mat'];
eegfile = [basename '.AccuSleep_EEG.mat'];
emgfile = [basename '.AccuSleep_EMG.mat'];

% load
load(emgfile, 'EMG')
load(eegfile, 'EEG')
load([basename, '.session.mat'])
load(labelsmanfile, 'labels')
gldstrd = labels;

% params
mousepath = fileparts(basepath);
SR = 512;
epochLen = 2.5;
minBoutLen = epochLen;
nstates = 4; 
labelnames = {'REM', 'WAKE', 'NREM', 'DROWSY'};

% separate data to 2
EMG_1st = EMG(1 : length(EMG) / 2);
EEG_1st = EEG(1 : length(EEG) / 2);
labels_1st = labels(1 : length(labels) / 2);

if mod(length(EMG), 2) ~= 0
    EMG_2nd = EMG(length(EMG) / 2 : length(EMG));
    EEG_2nd = EEG(length(EEG) / 2 : length(EEG));
else
    EMG_2nd = EMG(length(EMG) / 2 + 1 : length(EMG));
    EEG_2nd = EEG(length(EEG) / 2 + 1 : length(EEG));
end
if mod(length(labels), 2) ~= 0
    labels_2nd = labels(length(labels) / 2 : length(labels));
else
    labels_2nd = labels(length(labels) / 2 + 1 : length(labels));
end
if length(EMG) - length(EMG_1st) - length(EMG_2nd) ~= 0 
    warning('check data separation')
end
if length(labels) - length(labels_1st) - length(labels_2nd) ~= 0 
    warning('check data separation')
end

% visualize data
% AccuSleep_viewer(EEG, EMG, SR, epochLen, labels, [])
% AccuSleep_viewer(EEG_2nd, EMG_2nd, SR, epochLen, labels_2nd, [])
% AccuSleep_viewer(EEG_1st, EMG_1st, SR, epochLen, labels_1st, [])

% save vars
netpath = 'D:\Code\slutskycode\extracellular in vivo\lfp\detectStates\AccuSleep';
EEG = EEG_1st;
EMG = EMG_1st;
labels = labels_1st;
fileList{1, 1} = fullfile(netpath, [basename, '.AccuSleep_EEG1st.mat']);
fileList{1, 2} = fullfile(netpath, [basename, '.AccuSleep_EMG1st.mat']);
fileList{1, 3} = fullfile(netpath, [basename, '.AccuSleep_labels1st.mat']);
save(fileList{1, 1}, 'EEG')
save(fileList{1, 2}, 'EMG')
save(fileList{1, 3}, 'labels')

% train network
[net] = AccuSleep_train_nStates(fileList, SR, epochLen, 13, netpath);
netfile = 'D:\Code\AccuSleep\trainedNetworks\trainedNetwork2,5secEpochs.mat';
netfile = 'D:\Code\slutskycode\extracellular in vivo\lfp\detectStates\AccuSleep\4states_2,5s_6hrLabels_net.mat';
save('', 'net')         % !!! careful not to overwrite!!!
load(netfile, 'net')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% take data after training
EMG = EMG_2nd;
EEG = EEG_2nd;
gldstrd = labels_2nd;
nepochs = length(gldstrd);

% length of gldstrd labels used for creating the calibration data
calLen = nepochs;

% initialize
labelsOutput = zeros(length(calLen), nepochs);
labelsCal = ones(length(calLen), nepochs) * 4;
  
% create calibration 
calibrationData = createCalibrationData(standardizeSR(EEG, SR, 128),...
    standardizeSR(EMG, SR, 128), gldstrd, 128, epochLen);

% classify
labelsOutput = AccuSleep_classify(EEG, EMG, net, SR, epochLen,...
    calibrationData, minBoutLen);

% visualize data
AccuSleep_viewer(EEG, EMG, SR, epochLen, labelsOutput, [])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inspect results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
statesOutput = zeros(nstates, nstates);
edges = [1 : nstates + 1];
stateEpochs = histcounts(gldstrd);
for jj = 1 : nstates
    statesOutput(jj, :) = histcounts(labelsOutput(gldstrd == jj),...
        edges) / stateEpochs(jj) * 100;
end


AccuSimResults.calLen = calLen;
AccuSimResults.labelsCal = labelsCal;
AccuSimResults.labelsOutput = labelsOutput;
AccuSimResults.statesOutput = labelsOutput;
AccuSimResults.caldata = calibrationData;
AccuSimResults.labelnames = labelnames;
save(fullfile(mousepath, 'AccuSimResults.mat'), 'AccuSimResults')




