
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data base
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% hippocampal tetrodes ----------------------------------------------------
% acc
basepaths_acc = [{'G:\Data\lh99\as\lh99_211219_085802'},...     % local acsf. 120 min man
    {'G:\Data\lh99\as\lh99_211220_091903'},...                  % local ketamine. 150 min man
    {'G:\Data\lh93\lh93_210811_102035'},...                     % local ketamine. 150 min man
    {'D:\Data\lh93\lh93_210810_100810'},...                     % local ketamine. 120 min man
    {'D:\Data\lh93\lh93_210819_104212'},...                     % local saline. 70 min man
    {'D:\Data\lh93\lh93_210819_221608'},...                     % night. 70 min man
    {'I:\lh87\lh87_210523_100607'},...                          % op saline. 360 min man
    ];

test_acc = [{'G:\Data\lh98\lh98_211218_090630'},...             % local acsf. * min man
    {'I:\lh87\lh87_210522_095319'},...                          % op saline. 330 min man
    {'F:\Data\lh99\lh99_211218_090630'},...                     % local acsf. 120 min man
    ];

test_emg = [{'G:\Data\lh95\as\lh95_210825_080400'},...          % local ketamine. 120 min man
    {'G:\Data\lh99\as\lh99_211220_091903'},...                  % local ketamin. man done on acc
    {'F:\Data\Processed\lh96\lh96_211206_070400'},...           % ketamine i.p. 90 min man from 1st half, 30 min from 2nd half
    ];

% emg
basepaths_emg = [{'F:\Data\Processed\lh96\lh96_211201_070100'},...  % local acsf. 180 min man
    {'F:\Data\Processed\lh96\lh96_211205_072000'},...           % local mk801. 120 min man
    {'G:\Data\lh95\as\lh95_210824_083300'},...                  % local acsf. 200 min man                  
    {'I:\lh86\lh86_210302_183000'},...                          % op ketamine. 60 min man
    {'I:\lh86\lh86_210304_180800'},...                          % washout. 60 min man
    {'F:\Data\Processed\lh96\lh96_211207_071500'},...           % local baclofen. 60 min man
    ];

% frontal eeg -------------------------------------------------------------
% emg
basepaths = [{'I:\lh89\lh89_210510_085700'},...                 % op ketamine. 60 min man
    ];
% 'lh86'
% 'lh99'
% 'lh98'

% file list
fileList = as_fileLists(basepaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% config file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% general
fs = 1250;
epochLen = 1;
minBoutLen = epochLen;
nstates = 6; 

% load configuration file
[cfg_colors, cfg_names, cfg_weights, configfile] = as_loadConfig([]);

% weights
gldstrd = labels;
idx = gldstrd ~= 8;
weights = histcounts(gldstrd(idx)) / length(gldstrd(idx));
weights = round(weights * 100) / 100;       % round to two decimals
weights(4) = weights(4) + 1 - sum(weights); % add remainder to NREM
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
% state names
cfg_names = {'WAKE'; 'QWAKE'; 'LSLEEP'; 'NREM'; 'REM'; 'N/REM'; 'BIN'};

% save
save(configfile, 'cfg_colors', 'cfg_names', 'cfg_weights')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% train
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cntxEpochs = 63;
[net, netInfo] = AccuSleep_train(fileList, fs, epochLen, cntxEpochs);
netInfo.epochLen = epochLen;
netInfo.cntxEpochs = cntxEpochs;
netpath = 'D:\Code\slutsky_ECInVivo\lfp\SleepStates\AccuSleep\trainedNetworks';
netname = ['net_',  datestr(datetime, 'yymmdd_HHMMss')]; 
save(fullfile(netpath, netname), 'net', 'netInfo')      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test net (on remaining part of the data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% create calibration (based on entire dataset, doesn't matter)
calibrationData = createCalibrationData(standardizeSR(EEG_orig, fs, 128),...
    standardizeSR(EMG_orig, fs, 128), gldstrd, 128, epochLen);

% classify
[labels_net, netScores] = AccuSleep_classify(standardizeSR(EEG2, fs, 128),...
    standardizeSR(EMG2, fs, 128), net, 128, epochLen, calibrationData, minBoutLen);
    
% manually inspect model output and gldstrd
AccuSleep_viewer(EEG2, EMG2, fs, epochLen, labels_net, [])
AccuSleep_viewer(EEG2, EMG2, fs, epochLen, labels2, [])

% check precision
[netPrecision, netRecall] = as_cm(labels2, labels_net);


