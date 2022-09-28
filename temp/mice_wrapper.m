   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lh114
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% user input
<<<<<<< Updated upstream
mname           = 'lh111';
basepath        = 'I:\Data\lh111\2022-08-28_10-43-58';
=======
mname           = 'lh114';
basepath        = 'J:\Data\lh114\2022-09-26_10-00-36';
>>>>>>> Stashed changes
exp             = [1, 2];

% dat from oe 
mapch = 1 : 20;
rmvch = [5, 6];
rec = cell(max(exp), 1);
datInfo = preprocOE('basepath', basepath, 'exp', exp, 'rec', rec,...
    'rmvch', rmvch, 'mapch', mapch,...
    'nchans', length(mapch), 'fsIn', 20000);

% go to new folder
[mousepath, baseTime] = fileparts(basepath);
cd(mousepath)
dn = datenum(baseTime, 'yyyy-MM-dd');
recData = datestr(dn, 'yyMMdd');
fnames = dir(mousepath);
fidx = contains({fnames.name}, recData);
cd(fnames(fidx).name)

% session params
session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'forceDef', true, 'forceL', true, 'saveVar', true);      
basepath = session.general.basePath;
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;
[~, basename] = fileparts(basepath);

% clip bad parts
clear datfile
datfile = {fullfile(basepath, [basename, '.lfp'])};
clip = cell(1, length(datfile));
clip{1} = [seconds(minutes(575)), seconds(minutes(595));...
    seconds(minutes(517)), seconds(minutes(545))] * fs;
datInfo = preprocDat('orig_files', datfile, 'mapch', 1 : nchans,...
    'rmvch', [], 'nchans', nchans, 'saveVar', true,...
    'chunksize', 1e7, 'precision', 'int16', 'clip', clip);

% spike detection from temp_wh
[spktimes, ~] = spktimesWh('basepath', basepath, 'fs', fs, 'nchans', nchans,...
    'spkgrp', spkgrp, 'saveVar', false, 'saveWh', true,...
    'graphics', false, 'force', true, 'winWh', [0 Inf]);

% spike rate per tetrode
for igrp = 1 : length(spkgrp)
    spktimes{igrp} = spktimes{igrp} / fs;
end
sr = firingRate(spktimes, 'basepath', basepath,...
    'graphics', true, 'binsize', 60, 'saveVar', 'sr', 'smet', 'none',...
    'winBL', [0 Inf]);

% create ns files 
dur = [];
t = [];
spktimes2ns('basepath', basepath, 'fs', fs,...
    'nchans', nchans, 'spkgrp', spkgrp, 'mkClu', true,...
    'dur', dur, 't', t, 'grps', [1 : length(spkgrp)],...
    'spkFile', 'temp_wh');

% create lfp file
LFPfromDat('basepath', basepath, 'cf', 450, 'chunksize', 5e6,...
    'nchans', nchans, 'fsOut', 1250, 'fsIn', fs)    

% create emg signal from accelerometer data
acc = EMGfromACC('basepath', basepath, 'fname', [basename, '.lfp'],...
    'nchans', nchans, 'ch', nchans - 2 : nchans, 'saveVar', true, 'fsIn', 1250,...
    'graphics', false, 'force', true);

% sleep sig
sSig = as_prepSig([basename, '.lfp'], acc.mag,...
    'eegCh', [7 : 10], 'emgCh', [], 'saveVar', true, 'emgNchans', [],...
    'eegNchans', nchans, 'inspectSig', false, 'forceLoad', true,...
    'eegFs', 1250, 'emgFs', 1250, 'eegCf', [], 'emgCf', [10 450], 'fs', 1250);
labelsmanfile = [basename, '.sleep_labelsMan.mat'];
AccuSleep_viewer(sSig, [], labelsmanfile)

% calc spec
spec = calc_spec('sig', [], 'fs', 1250, 'graphics', true,...
    'saveVar', true, 'padfft', -1, 'winstep', 5,...
    'ftarget', [], 'ch', {[7 : 10]},...
    'force', true, 'logfreq', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lh112
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% user input
mname           = 'lh112';
<<<<<<< Updated upstream
basepath        = 'D:\Data\lh112\2022-09-03_11-23-20';
=======
basepath        = 'J:\Data\lh112\2022-09-03_11-23-20';
>>>>>>> Stashed changes
exp             = [1];

% dat from oe 
mapch = 1 : 20;
rmvch = [];
rec = cell(max(exp), 1);
datInfo = preprocOE('basepath', basepath, 'exp', exp, 'rec', rec,...
    'rmvch', rmvch, 'mapch', mapch,...
    'nchans', length(mapch), 'fsIn', 20000);

% go to new folder
[mousepath, baseTime] = fileparts(basepath);
cd(mousepath)
dn = datenum(baseTime, 'yyyy-MM-dd');
recData = datestr(dn, 'yyMMdd');
fnames = dir(mousepath);
fidx = contains({fnames.name}, recData);
cd(fnames(fidx).name)

% session params
session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'forceDef', true, 'forceL', true, 'saveVar', true);      
basepath = session.general.basePath;
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;
[~, basename] = fileparts(basepath);

% clip bad parts
clear datfile
datfile = {fullfile(basepath, [basename, '.dat'])};
clip = cell(1, length(datfile));
clip{1} = [seconds(minutes(755)), Inf] * fs;
datInfo = preprocDat('orig_files', datfile, 'mapch', 1 : nchans,...
    'rmvch', [], 'nchans', nchans, 'saveVar', true,...
    'chunksize', 1e7, 'precision', 'int16', 'clip', clip);

% spike detection from temp_wh
[spktimes, ~] = spktimesWh('basepath', basepath, 'fs', fs, 'nchans', nchans,...
    'spkgrp', spkgrp, 'saveVar', true, 'saveWh', true,...
    'graphics', false, 'force', true, 'winWh', [0 Inf]);

% spike rate per tetrode
for igrp = 1 : length(spkgrp)
    spktimes{igrp} = spktimes{igrp} / fs;
end
sr = firingRate(spktimes, 'basepath', basepath,...
    'graphics', true, 'binsize', 60, 'saveVar', 'sr', 'smet', 'none',...
    'winBL', [0 Inf]);

% create ns files 
dur = [];
t = [];
spktimes2ns('basepath', basepath, 'fs', fs,...
    'nchans', nchans, 'spkgrp', spkgrp, 'mkClu', true,...
    'dur', dur, 't', t, 'grps', [1 : length(spkgrp)],...
    'spkFile', 'temp_wh');

% create lfp file
LFPfromDat('basepath', basepath, 'cf', 450, 'chunksize', 5e6,...
    'nchans', nchans, 'fsOut', 1250, 'fsIn', fs)    

% create emg signal from accelerometer data
acc = EMGfromACC('basepath', basepath, 'fname', [basename, '.lfp'],...
    'nchans', nchans, 'ch', nchans - 2 : nchans, 'saveVar', true, 'fsIn', 1250,...
    'graphics', false, 'force', true);

% sleep sig
sSig = as_prepSig([basename, '.lfp'], acc.mag,...
    'eegCh', [7], 'emgCh', [], 'saveVar', true, 'emgNchans', [],...
    'eegNchans', nchans, 'inspectSig', false, 'forceLoad', true,...
    'eegFs', 1250, 'emgFs', 1250, 'eegCf', [], 'emgCf', [10 450], 'fs', 1250);
labelsmanfile = [basename, '.sleep_labelsMan.mat'];
AccuSleep_viewer(sSig, [], labelsmanfile)

% calc spec
spec = calc_spec('sig', [], 'fs', 1250, 'graphics', true,...
    'saveVar', true, 'padfft', -1, 'winstep', 5,...
    'ftarget', [], 'ch', {[7]},...
    'force', true, 'logfreq', true);


