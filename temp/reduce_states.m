% reduce number of states 

% manual paths for 6 states
basepaths{1} = 'F:\Data\Colleagues\HB\sleepStates\WT1\baseline_210419_0822';
basepaths{2} = 'F:\Data\Colleagues\HB\sleepStates\Old_wt2\wt2_200228_0700';
basepaths{3} = 'F:\Data\Colleagues\HB\sleepStates\Old_wt1\WT1_200327_1047';
basepaths{4} = 'F:\Data\Colleagues\HB\sleepStates\WT3\baseline_210622_2015';
netfile = 'D:\Code\slutskycode\extracellular in vivo\lfp\SleepStates\AccuSleep\trainedNetworks\net_211107_141044.mat';

% manual paths for 4 states
basepaths{1} = 'F:\Data\Colleagues\HB\sleepStates\WT1\4states';
basepaths{2} = 'F:\Data\Colleagues\HB\sleepStates\Old_wt2\4states';
basepaths{3} = 'F:\Data\Colleagues\HB\sleepStates\Old_wt1\4states';
basepaths{4} = 'F:\Data\Colleagues\HB\sleepStates\WT3\4states';
netfile = 'D:\Code\slutskycode\extracellular in vivo\lfp\SleepStates\AccuSleep\trainedNetworks\net_211104_171343.mat';


% selected states, the rest will be overwritten
statesReplace = [3, 4; 6, 4; 4, 3; 5, 4; 7, 5; 8, 6];

for ipath = 1 : length(basepaths)
    
    % load
    cd(basepaths{ipath})
    eegfile = dir('*EEG*');
    load(eegfile(1).name)
    emgfile = dir('*EMG*');
    load(emgfile(1).name)
    labelsfile = dir('*labelsMan*');
    load(labelsfile(1).name)
    
    % reduce to nsatates
    for istate = 1 : size(statesReplace, 1)
        labels(labels == statesReplace(istate, 1)) =...
            statesReplace(istate, 2);
    end
    
    % save
    newpath{ipath} = fullfile(basepaths{ipath}, '4states');
    mkdir(newpath{ipath})
    cd(newpath{ipath})
    [~, basename] = fileparts(newpath{ipath});
    
    labelsmanfile = [basename, '.AccuSleep_labelsMan.mat'];
    eegfile = [basename '.AccuSleep_EEG.mat'];
    emgfile = [basename '.AccuSleep_EMG.mat'];
        
    save(labelsmanfile, 'labels')
    save(eegfile, 'EEG')
    save(emgfile, 'EMG')
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% config file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load configuration file
[cfg_colors, cfg_names, cfg_weights, configfile] = as_loadConfig([]);
config_bkup = [fileparts(configfile), '\AS_config_6states'];
save(config_bkup, 'cfg_colors', 'cfg_names', 'cfg_weights')

% weights
clear cfg_weigths
gldstrd = labels;
idx = gldstrd ~= 6;
weights = histcounts(gldstrd(idx)) / length(gldstrd(idx));
weights = round(weights * 100) / 100;       % round to two decimals
weights(3) = weights(3) + 1 - sum(weights); % add remainder to NREM
cfg_weights = weights;

% colors
clear cfg_colors
cfg_colors{1} = [240 110 110] / 255;
cfg_colors{2} = [240 170 125] / 255;
cfg_colors{3} = [110 180 200] / 255;
cfg_colors{4} = [170 100 170] / 255;
cfg_colors{5} = [200 200 200] / 255;
cfg_colors = cfg_colors(:);

% state names
cfg_names = {'WAKE'; 'QWAKE'; 'NREM'; 'REM'; 'BIN'};

% save new config
save(configfile, 'cfg_colors', 'cfg_names', 'cfg_weights')

% validate new config
AccuSleep_viewer(EEG, EMG,  1250, 1, labels, [])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% train
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% arange filelists
fileList = as_fileLists(basepaths);

% general
fs = 1250;
boutLen = 1;
minBoutLen = boutLen;

[net, netInfo] = AccuSleep_train(fileList, fs, boutLen, 63);
netpath = 'D:\Code\slutskycode\extracellular in vivo\lfp\SleepStates\AccuSleep\trainedNetworks';
netname = ['net_',  datestr(datetime, 'yymmdd_HHMMss')]; 
save(fullfile(netpath, netname), 'net', 'netInfo')      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path to test data
testpath = 'F:\Data\Colleagues\HB\sleepStates\WT3\test';
cd(testpath)

% load
eegfile = dir('*EEG*');
load(eegfile(1).name)
emgfile = dir('*EMG*');
load(emgfile(1).name)
labelsfile = dir('*labelsMan*');
load(labelsfile(1).name)

% netfile = [];
ss = as_wrapper(EEG, EMG, [], 'basepath', pwd, 'calfile', [],...
    'viaGui', false, 'forceCalibrate', true, 'inspectLabels', false,...
    'saveVar', true, 'forceAnalyze', true, 'fs', 1250, 'netfile', netfile,...
    'graphics', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examine netScores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gldstrd = labels;

[netPrecision, netRecall] = as_cm(gldstrd, ss.labels_net, ss.netScores,...
    'graphics', true);


setMatlabGraphics(false)
newLabels = ss.labels_net;
nstates = size(ss.netScores, 2);
thr = 0 : 0.1 : 1;

for ithr = 1 : length(thr)      
    for istate = 1 : nstates
        
        stateIdx = ss.labels_net == istate;
        scoreIdx = ss.netScores(:, istate) < thr(ithr) & stateIdx;        
        lostData(ithr, istate) = sum(scoreIdx) / sum(stateIdx) * 100;
        newLabels(scoreIdx) = nstates + 2;                   
    end

    [netPrecision(ithr, :), netRecall(ithr, :)] =...
        as_cm(gldstrd, newLabels, 'graphics', false);    
end

% graphics ----------------------------------------------------------------

fh = figure;
for istate = 1 : nstates
    
    subplot(2, nstates / 2, istate)
    plot(thr, netPrecision(:, istate) * 100, 'k', 'LineWidth', 2)
    hold on
    plot(thr, netRecall(:, istate) * 100, 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
    xlabel('Threshold')
    ylabel('Performace [%]')
    ylim([50 100])
    
    yyaxis right
    ph = plot(thr, lostData(:, istate), 'LineWidth', 2);
    ph.Color = cfg_colors{istate};
    ylabel('Data Lost [%]')
    set(gca, 'ycolor', cfg_colors{istate})
    ylim([0 100])
    
    set(gca, 'box', 'off', 'TickLength', [0 0])
    title(cfg_names{istate})
    if istate == 1
        legend({'Precision', 'Recall', 'DataLost'})
    end
end



