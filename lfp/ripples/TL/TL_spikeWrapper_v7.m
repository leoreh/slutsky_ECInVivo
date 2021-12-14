function TL_spikeWrapper_v7(basepath , chooseSet , runDat , runLFP , runLFPnoFilt , runWhite , runEMG , runSleep , runSU)

% FUNCTION:
% Wrapper code to convert raw tdt files to dat, lfp, spiketimes, and single
% units. _v7 allows user to select which files to produce (see inputs).

% Folder/File structure: 
% (1) Should have folder titled 'Dates', within which
% are folders named for the start date of each recording YYMMDD (eg
% '210825','210826' ..). 
% (2) Within each dated folder must be an excel spreadsheet with experiment 
% information. Each excel spreadsheet has a tab  for each subject.
% (3) Code will create a folder for each subject found in a 'Subjects' 
% folder adjacent to the 'Dates' folder and save the date under the 
% appropriate subject and date folder..

% eg TL_spikeWrapper('E:\Tomer Spikes\Dates\210610') run on subjects 'MCU1'
% and 'MCU2' will save data into ->
% 'E:\Tomer Spikes\Subjects\MCU1\210610' and
% 'E:\Tomer Spikes\Subjects\MCU2\210610'

% INPUTS:
% [basepath] : string path to recording date (contains folders for each
% block on that date and an excel spreadsheet with recording info)
% [chooseSet] : set number to run (same set for all subjects as of _v7
% [runDat; runLFP; runLFPnofilt; runWhite; runEMG; runSleep; runSU] : 0 to
% run for all subjects, 1 for Raw1, 2 for Raw2
%% Read excel spreadsheet
% There must be an excel spreadsheet in the basepath (..\Dates\YYMMDD)
% containing experiment information. This section will convert this to
% expInfo, a mat file that will be read in all later functions. tdt2dat
% code will further edit expInfo with timing parameters extracted from tdt.

% correct basepath name if necessary
if ~strcmp(basepath(end) , filesep)
    basepath = [basepath , filesep];
end

% Load excel with experiment information
% [subjects] : cell array of subject names in excel file
% [expath] : string path to excel file
temp = dir(basepath);
fnames = {temp(~[temp.isdir]).name}';
clear temp;
temp = find(cell2mat(cellfun(@(z) ~isempty(strfind(z,'xls')),fnames,'uniformoutput',false)));
expath = [basepath  , fnames{temp}];
[~ , subjects] = xlsfinfo(expath)
clear temp fnames;

% Loop thru each subject
for s = 1 : length(subjects)
    
    % Load excel tab for current subject
    [aa , bb , excel] = xlsread(expath , subjects{s});
    headers = [excel(1,:)];
    
% [expInfo] : structure containing experiment info. This section will just
% rip the info directly off the excel spreadsheet. The expInfo will be
% updated in tdt2dat with length of recording and true start time
    for r = 2:size(excel , 1)
        for h = 1 : length(headers)
            val = [];
            og = excel{r,h};
            if isnumeric(og) % is it nan
                val = og;
            else
                s2n = str2num(og);
                if ~isempty(s2n) % is it a number
                    val = s2n;
                else
                    sf = strfind(og , ',');
                    if ~isempty(sf) % is it strings separated by commas
                        sf = [0 , sf];
                        for ssf = 1 : length(sf)
                            if ssf < length(sf)
                                val{ssf} = og(sf(ssf)+1:sf(ssf+1)-1);
                            else
                                val{ssf} = og(sf(ssf)+1:end);
                            end
                        end
                    else
                        val = og;
                    end
                end
                
            end
            expInfo.(headers{h}){r-1,1} = val;
            clear og val;
        end
    end
    clear r h;
    
    % 1 field in expInfo is "set" which says which blocks should be sorted
    % together (as a set)
    % set up folders for saving data into

    allSets = [expInfo.set{:}]; % set value for each block
    [uSets i] = unique(allSets);
    fs = strfind(basepath , filesep);
    dates = expInfo.date(i);
    for d = 1 : length(dates)
        saveFol = [basepath(1:fs(end-2)) 'Subjects' , filesep , subjects{s} , filesep , num2str(dates{d})];
        expInfo.saveFol{d} = saveFol;
        if ~isfolder(saveFol)
            mkdir(saveFol);
        end
    end
    clear i fs dates d saveFol;
    
    %% Convert TDT 2 Dat for the current subject
    
    % Now loop thru each set of blocks and convert from tdt to dat
    if isempty(chooseSet)
        setRange = 1 : length(uSets);
    else
        setRange = find(uSets == chooseSet);
    end
    
    for us = setRange
        blockI = find(allSets == uSets(us)); % find which blocks correspond to this set
%         us=2;
        if ~isempty(runDat)
            if strcmp(expInfo.store{blockI(1)} , 'Raw1') & runDat == 1 | ...
                    strcmp(expInfo.store{blockI(1)} , 'Raw2') & runDat == 2 | ...
                    runDat == 0
                expInfo = TL_tdt2dat_v3(expInfo , uSets(us));
                % user needs to create XML file to continue:
                % 1. open dat in neuroscope and close it (XML file will be made)
                % 2. open XML file, go to spike groups and select channels per
                % tetrode group (0 indexed), set samples per waveform to 40, and
                % set # features to 3
                display(['Completed dat file for ' expInfo.subject{blockI(1)} ' ' num2str(expInfo.date{blockI(1)})]);
                display(['Remove channel ' num2str(expInfo.rmvch{blockI(1)})]);
                display('After making XML file and setting spike groups in ND Manager, press any key to continue');
                pause;
            end
        else
            load([expInfo.saveFol{blockI(1)} , filesep expInfo.store{blockI(1)} , '_expInfo.mat']);
        end
        
        
        if ~isempty(runLFP)
            if strcmp(expInfo.store{blockI(1)} , 'Raw1') & runLFP == 1 | ...
                    strcmp(expInfo.store{blockI(1)} , 'Raw2') & runLFP == 2 | ...
                    runLFP == 0 ,
                
                % Make LFP file sampled at 1250 Hz
                % TL EDIT INPUT TO NCHANS AND FSIN TO USE FROM EXPINFO
                TL_LFPfromDat(structfun(@(z) z{blockI(1)} , expInfo , 'uniformoutput' , false) , 'basepath' , expInfo.saveFol{blockI(1)} , 'cf' , 450 , 'chunksize' , 5e6 ,...
                    'nchans', length(expInfo.channels{blockI(1)}), 'fsOut', 1250,...
                    'fsIn', 24414.06)
            end
        end
        
        if ~isempty(runLFPnoFilt)
            if strcmp(expInfo.store{blockI(1)} , 'Raw1') & runLFPnoFilt == 1 | ...
                    strcmp(expInfo.store{blockI(1)} , 'Raw2') & runLFPnoFilt == 2 | ...
                    runLFPnoFilt == 0 ,
                
                % Make LFP file no filter less downsample (used to
                % visualize raw ripples (cant see spikes with 250 cf)
                TL_LFPfromDat(structfun(@(z) z{blockI(1)} , expInfo , 'uniformoutput' , false) , 'basepath' , expInfo.saveFol{blockI(1)} , 'cf' , 8000 , 'chunksize' , 5e6 ,...
                    'nchans', length(expInfo.channels{blockI(1)}), 'fsOut', 12500,...
                    'fsIn', 24414.06)
            end
        end
        
        if ~isempty(runWhite)
            if strcmp(expInfo.store{blockI(1)} , 'Raw1') & runWhite == 1 | ...
                    strcmp(expInfo.store{blockI(1)} , 'Raw2') & runWhite == 2 | ...
                    runWhite == 0 ,
                
                fs = 24414.06;
                
                % convert channel mapping to be consecutive with rmvch
                % removed
                rmvch = expInfo.rmvch{blockI(1)};
                fullSpkgrp = {[1:4],[5:8],[9:12],[13:16]};
                temp = cell2mat(cellfun(@(z) 4 - sum(ismember(z,rmvch)) , fullSpkgrp , 'uniformoutput' , false));
                sm = 0;
                for t = 1 : length(temp)
                    spkgrp{t} = [1:temp(t)] + sm;
                    sm = sm + temp(t);
                end
                clear temp;
                nchans = length(expInfo.channels{blockI(1)});% +  length(expInfo.rmvch{blockI(1)});
                
                [spktimes, ~] = TL_spktimesWh('basepath', expInfo.saveFol{blockI(1)}, 'fs', fs, 'nchans', nchans,...
                    'spkgrp', spkgrp, 'saveVar', true, 'saveWh', true,...
                    'graphics', false, 'force', true);
                
            end
        end
        
        
        if ~isempty(runEMG)
            if strcmp(expInfo.store{blockI(1)} , 'Raw1') & runEMG == 1 | ...
                    strcmp(expInfo.store{blockI(1)} , 'Raw2') & runEMG == 2 | ...
                    runEMG == 0 ,
                
                %                 Extract EMG
                emgInfo = expInfo;
                emgInfo.store{(uSets(us))} = ['EMG' , expInfo.store{uSets(us)}(end)];
                emgInfo.rmvch{(uSets(us))} = [];
                emgInfo.channels{(uSets(us))} = 1:4;
                emgInfo.fs = [];
                
                emgInfo = TL_tdt2dat_v3(emgInfo , uSets(us));
            end
            %                         TL_tdt2dat([emgInfo(blockI)]);
            
            
        end
        
        
        
        if ~isempty(runSleep)
            
            if strcmp(expInfo.store{blockI(1)} , 'Raw1') & runSleep == 1 | ...
                    strcmp(expInfo.store{blockI(1)} , 'Raw2') & runSleep == 2 | ...
                    runSleep == 0 ,
                
                temp = dir([expInfo.saveFol{blockI(1)} , filesep , '*EMG' , '*.dat']);
                emgFname = [temp.folder , filesep , temp.name];
                
                temp = dir([expInfo.saveFol{blockI(1)} , filesep , '*EMG' , '*expInfo.mat']);
                emgInfo = load([temp.folder , filesep , temp.name]);
                emgInfo = emgInfo.expInfo;
%                 emgInfo = expInfo;
                
                lfpFname = [expInfo.saveFol{blockI(1)} , filesep , expInfo.store{blockI(1)} , '_' , 'Set' , num2str(expInfo.set{blockI(1)}) , '_' , num2str(expInfo.date{blockI(1)}) , '.lfp'];
                % TL emgCh should be entered into column of excel spreadsheet so need
                % to select it during the recording. Otherwise can use a pause and wait
                % for user input
                
                [EMG , EEG , sigInfo] = TL_as_prepSig(lfpFname , emgFname , expInfo , uSets(us) , ...
                    'eegCh' , [6] , 'emgCh' , 1 , 'saveVar' , true , 'emgNchans', 4 , 'eegNchans' ,  length(expInfo.channels{blockI(1)}) , ...
                    'inspectSig' , true , 'forceLoad' , true, 'eegFs' , 1250 , 'emgFs' ,emgInfo.fs{blockI(1)}) ;
                %emgInfo.fs{blockI(1)} ,
            end
        end
        
                if ~isempty(runSU)
            
            if strcmp(expInfo.store{blockI(1)} , 'Raw1') & runSU == 1 | ...
                    strcmp(expInfo.store{blockI(1)} , 'Raw2') & runSU == 2 | ...
                    runSU == 0 ,
             
                % convert channel mapping to be consecutive with rmvch
                % removed
                rmvch = expInfo.rmvch{blockI(1)};
                fullSpkgrp = {[1:4],[5:8],[9:12],[13:16]};
                temp = cell2mat(cellfun(@(z) 4 - sum(ismember(z,rmvch)) , fullSpkgrp , 'uniformoutput' , false));
                sm = 0;
                for t = 1 : length(temp)
                    spkgrp{t} = [1:temp(t)] + sm;
                    sm = sm + temp(t);
                end
                clear temp;
                
                TL_spktimes2ns(expInfo , blockI , 'basepath', expInfo.saveFol{blockI(1)}, 'fs', expInfo.fs{blockI(1)},...
    'nchans', length(expInfo.channels{blockI(1)}), 'spkgrp', spkgrp, 'mkClu', true,...
    'dur', [], 't', [], 'psamp', [], 'grps', [1 : length(spkgrp)],...
    'spkFile', 'temp_wh');

                
            end
        end
        
        
    end
    
    
    
    
end

end
%
%
%
%
%         % Run this function to make 'session' which includes experiment
%         % information, but not really sure (legacy) .. it is needed for
%         % some functions below I believe...
% %         try
% %         session = TL_CE_sessionTemplate(saveFol, expInfo , 'viaGUI', false,...
% %             'force', true, 'saveVar', true);
% %         catch
% %             display('Try making XML and spike groups again..you screwed up. Then press any key to continue');
% %             pause;
% %             session = TL_CE_sessionTemplate(saveFol, expInfo , 'viaGUI', false,...
% %             'force', true, 'saveVar', true);
% %         end
% %         basepath = session.general.basePath;
% %         nchans = session.extracellular.nChannels;
% %         fs = session.extracellular.sr;
% %         spkgrp = session.extracellular.spikeGroups.channels;
% %         [~, basename] = fileparts(basepath);
%
%% Extract spike times with TL_spktimesWh

%
%         % fix spkgrp to match spktimes2ns code (accomodate for removed
%         % channels..)
%         tot = 0;
%         for s = 1 : length(spkgrp)
%             sg{s} = [1:length(spkgrp{s})] + tot;
%             tot = sg{s}(end);
%         end
%         dur = [];  % TL STILL DO NOT CONNECT THIS WITH EXPINFO
%         t = [];  % TL STILL DO NOT CONNECT THIS WITH EXPINFO
%         TL_spktimes2ns(expInfo , uSets(us) , 'basepath' , saveFol , ...
%             'fs' , fs , 'nchans' , nchans , 'spkgrp', sg , 'mkClu' , true ,...
%             'dur', dur , 't', t , 'psamp', [] , 'grps' , [1 : length(spkgrp)] ,...
%             'spkFile' , 'temp_wh');
%
%         % Extract EMG
%         emgInfo = expInfo;
%         for e = 1 : length(expInfo{us})
%             emgInfo{e}.store = ['EMG' , expInfo{e}.store(end)];
%             emgInfo{e}.rmvch = [];
%             emgInfo{e}.channels = 1:4;
%         end
%         TL_tdt2dat([emgInfo(blockI)]);
%     end
%     end




















%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % session info (cell explorer foramt)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% session = CE_sessionTemplate(pwd, 'viaGUI', false,...
%     'force', true, 'saveVar', true);
% basepath = session.general.basePath;
% nchans = session.extracellular.nChannels;
% fs = session.extracellular.sr;
% spkgrp = session.extracellular.spikeGroups.channels;
% [~, basename] = fileparts(basepath);
%
% %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % spike sorting
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % spike detection from temp_wh
%
%
% % spike rate
% for ii = 1 : length(spkgrp)
%     spktimes{ii} = spktimes{ii} / fs;
% end
% sr = firingRate(spktimes, 'basepath', basepath,...
%     'graphics', false, 'saveFig', false,...
%     'binsize', 60, 'saveVar', 'sr', 'smet', 'none',...
%     'winBL', [0 Inf]);
%
% % clip ns files
% % dur = -420;
% % t = [];
% % nsClip('dur', dur, 't', t, 'bkup', true, 'grp', [3 : 4]);
%
% create ns files
% dur = [];
% t = [];
% spktimes2ns('basepath', basepath, 'fs', fs,...
%     'nchans', nchans, 'spkgrp', spkgrp, 'mkClu', true,...
%     'dur', dur, 't', t, 'psamp', [], 'grps', [1 : length(spkgrp)],...
%     'spkFile', 'temp_wh');
%
% % clean clusters after sorting
% cleanCluByFet('basepath', pwd, 'manCur', false, 'grp', [1 : 4])
%
% % cut spk from dat and realign
% fixSpkAndRes('grp', 3, 'dt', 0, 'stdFactor', 0);
%
% % cell explorer metrics
% cell_metrics = ProcessCellMetrics('session', session,...
%     'manualAdjustMonoSyn', false, 'summaryFigures', false,...
%     'debugMode', true, 'transferFilesFromClusterpath', false,...
%     'submitToDatabase', false, 'getWaveformsFromDat', true);
% cell_metrics = CellExplorer('basepath', basepath);
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % spikes
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % cluster validation
% mu = [];
% spikes = cluVal('spikes', spikes, 'basepath', basepath, 'saveVar', true,...
%     'saveFig', false, 'force', true, 'mu', mu, 'graphics', false,...
%     'vis', 'on', 'spkgrp', spkgrp);
%
% % firing rate
% binsize = 60;
% winBL = [1 Inf];
% fr = firingRate(spikes.times, 'basepath', basepath, 'graphics', false, 'saveFig', false,...
%     'binsize', binsize, 'saveVar', true, 'smet', 'MA', 'winBL', winBL);
%
% % CCG
% binSize = 0.001; dur = 0.12; % low res
% binSize = 0.0001; dur = 0.02; % high res
% [ccg, t] = CCG({xx.times{:}}, [], 'duration', dur, 'binSize', binSize);
% u = 20;
% plotCCG('ccg', ccg(:, u, u), 't', t, 'basepath', basepath,...
%     'saveFig', false, 'c', {'k'}, 'u', spikes.UID(u));
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % lfp
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % create lfp
% LFPfromDat('basepath', basepath, 'cf', 450, 'chunksize', 5e6,...
%     'nchans', nchans, 'fsOut', 1250,...
%     'fsIn', fs)
%
% % load lfp
% lfpInterval = [150 * 60, 180 * 60];
% lfp = getLFP('basepath', basepath, 'ch', [spkgrp{:}], 'chavg', {},...
%     'fs', 1250, 'interval', lfpInterval, 'extension', 'lfp',...
%     'savevar', true, 'forceL', true, 'basename', '');
%
% % remove 50 Hz from signal
% emgOrig = filterLFP(emgOrig, 'fs', 1250, 'stopband', [49 51],...
%     'dataOnly', true, 'saveVar', false, 'graphics', false);
%
% acc = EMGfromACC('basepath', basepath, 'fname', [basename, '.lfp'],...
%     'nchans', 20, 'ch', [17 : 19], 'saveVar', true, 'fsIn', 1250,...
%     'graphics', false);
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % sleep states
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % prep signal (alt 2)
% [EMG, EEG, sigInfo] = as_prepSig([basename, '.lfp'], acc.mag,...
%     'eegCh', [9 : 12], 'emgCh', [], 'saveVar', true, 'emgNchans', [],...
%     'inspectSig', true, 'forceLoad', true, 'eegFs', 1250, 'emgFs', 1250);
%
% % manually create labels
% labelsmanfile = [basename, '.AccuSleep_labelsMan.mat'];
% AccuSleep_viewer(EEG, EMG, 1250, 1, [], labelsmanfile)
%
% % classify with a network
% ss = as_wrapper(EEG, EMG, [], 'basepath', basepath, 'calfile', [],...
%     'viaGui', false, 'forceCalibrate', true, 'inspectLabels', true,...
%     'saveVar', true, 'forceAnalyze', true, 'fs', 1250);
%
% % inspect separation after manual scoring
% as_inspectSeparation(EEG, EMG, labels)
%
% % get confusion matrix between two labels
% [ss.netPrecision, ss.netRecall] = as_cm(labels1, labels2);
%
% % show only x hours of data
% x = 11;
% tidx = [1 : x * 60 * 60 * 1250];
% lidx = [1 : x * 60 * 60];
% AccuSleep_viewer(EEG(tidx), EMG(tidx), 1250, 1, labels(lidx), [])
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % handle dat  files
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %%% cat dat
% nchans = 17;
% fs = 1250;
% saveFol = mousepath;
% datFiles{1} = 'K:\Data\lh91\experiment11\recording1\continuous\Rhythm_FPGA-109.0\continuous.dat';
% datFiles{2} = 'K:\Data\lh91\experiment11\recording1\continuous\Rhythm_FPGA-109.0\28e1.dat';
% sigInfo = dir(datFiles{1});
% nsamps = floor(sigInfo.bytes / class2bytes('int16') / nchans);
% parts{1} = [nsamps - 2 * 60 * 60 * fs nsamps];
% parts{2} = [0 4 * 60 * 60 * fs];
%
% catDatMemmap('datFiles', datFiles, 'saveFol', saveFol, 'parts', parts,...
%     'nchans', nchans, 'saveVar', true)
%
%
% % preproc dat
% clip = [1, 864000000];
% datInfo = preprocDat('basepath', pwd,...
%     'fname', 'continuous.dat', 'mapch', 1 : 20,...
%     'rmvch', [3, 7, 13], 'nchans', 20, 'saveVar', false, 'clip', clip,...
%     'chunksize', 5e6, 'precision', 'int16', 'bkup', true);
%
%
