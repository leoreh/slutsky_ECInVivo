% accusleep_simulation. trains a network on half a data set and test
% accuracy on the other half

% step 1 - configuration
% step 2 - prepare data files (e.g. split 2 two)
% step 3 - train network on 1st half
% step 4 - create calibration on 2nd half using gldstrd
% step 5 - classify 2nd half
% step 6 - compare output with gldstrd
% step 7 - save output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% general
fs = 1250;
epochLen = 1;
minBoutLen = epochLen;
nstates = 6; 

% load configuration file
configfile = 'D:\Code\slutskycode\extracellular in vivo\lfp\SleepStates\AccuSleep\AS_config.mat';
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

% state names
cfg_names = {'WAKE'; 'QWAKE'; 'LSLEEP'; 'NREM'; 'REM'; 'N/REM'; 'BIN'};

% calc weights
idx = gldstrd ~= 8;
weights = histcounts(gldstrd(idx)) / length(gldstrd(idx));
weights = round(weights * 100) / 100;       % round to two decimals
weigths = [weights 0];
weights(4) = weights(4) + 1 - sum(weights); % add remainder to NREM
% weights = [0.32 0.16 0.06 0.38 0.08 0 0];   % overwrite
cfg_weights = weights;
cfg_weights = [cfg_weights 0];

% colors
cfg_colors{1} = [240 110 110] / 255;
cfg_colors{2} = [240 170 125] / 255;
cfg_colors{3} = [150 205 130] / 255;
cfg_colors{4} = [110 180 200] / 255;
cfg_colors{5} = [170 100 170] / 255;
cfg_colors{6} = [200 200 100] / 255;
cfg_colors{7} = [200 200 200] / 255;

% save
save(configfile, 'cfg_colors', 'cfg_names', 'cfg_weights')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepaths{1} = 'E:\Data\lh86\lh86_210227_190600';
cd(basepaths{1})
[~, basename] = fileparts(basepaths{1});

labelsmanfile = [basename, '.AccuSleep_labelsMan.mat'];
eegfile = [basename '.AccuSleep_EEG.mat'];
emgfile = [basename '.AccuSleep_EMG.mat'];

% load
load(emgfile, 'EMG')
load(eegfile, 'EEG')
load([basename, '.session.mat'])
load(labelsmanfile, 'labels')
gldstrd = labels;

% ALT 1: train network on entire data
fileList = as_fileLists(basepaths);


netpath = 'D:\Code\slutskycode\extracellular in vivo\lfp\SleepStates\AccuSleep\trainedNetworks';
netname = 'lh86_4hr'; 
[net] = AccuSleep_train(fileList, fs, epochLen, 13, 'E:\');
save(fullfile(netpath, netname), 'net')         % !!! careful not to overwrite!!!

% ALT 2: separate data to 2
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
% AccuSleep_viewer(EEG, EMG, fs, epochLen, labels, [])
AccuSleep_viewer(EEG_2nd, EMG_2nd, fs, epochLen, labels_2nd, [])
AccuSleep_viewer(EEG_1st, EMG_1st, fs, epochLen, labels_1st, [])

% save vars
netpath = 'D:\Data\ser2';
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
[net, trainInfo] = AccuSleep_train(fileList, fs, epochLen, 13, netpath);
save('D:\Data\ser2\6states_1s_6hrLabels_HBser2_net', 'net')         % !!! careful not to overwrite!!!
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
calibrationData = createCalibrationData(standardizefs(EEG, fs, 128),...
    standardizefs(EMG, fs, 128), gldstrd, 128, epochLen);

% classify
[labelsOutput, scores] = AccuSleep_classify(EEG, EMG, net, fs, epochLen,...
    calibrationData, minBoutLen);

% visualize data
AccuSleep_viewer(EEG, EMG, fs, epochLen, labelsOutput, [])
AccuSleep_viewer(EEG, EMG, fs, epochLen, gldstrd, [])

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inspect results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
statesOutput = zeros(nstates, nstates);
edges = [1 : nstates + 1];
% accuracy
stateEpochs = histcounts(gldstrd);
for istate = 1 : nstates
    statesOutput(istate, :) = histcounts(labelsOutput(gldstrd == istate),...
        edges) / stateEpochs(istate) * 100;
end

stateEpochs = histcounts(labelsClean, edges);
for istate = 1 : nstates
    statesOutput(istate, :) = histcounts(gldstrd(labelsClean == istate),...
        edges) / stateEpochs(istate) * 100;
end

AccuSimResults.calLen = calLen;
AccuSimResults.labelsCal = labelsCal;
AccuSimResults.labelsOutput = labelsOutput;
AccuSimResults.statesOutput = labelsOutput;
AccuSimResults.caldata = calibrationData;
AccuSimResults.labelnames = cfg_names;
save(fullfile(mousepath, 'AccuSimResults.mat'), 'AccuSimResults')




