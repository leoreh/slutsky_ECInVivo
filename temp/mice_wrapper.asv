    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lh106
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% user input
mname = 'lh106';
basepath = 'I:\Data\lh106\lh106_220525_091000';
blocks = [1 : 2];

% tank to dat (tetrodes)
store = 'Raw1';
mapch = 1 : 16;
rmvch = [];
chunksize = 300;
clip = cell(1, 1);
datInfo = tdt2dat('basepath', basepath, 'store', store, 'blocks',  blocks,...
    'chunksize', chunksize, 'mapch', mapch, 'rmvch', rmvch, 'clip', clip);

% move dat to session folder
[~, basename] = fileparts(basepath);
recname = strrep(basename, 'lh106', mname);
recpath = fullfile(basepath, recname);
mkdir(recpath)
fnames = dir(['*' basename '*']);
for ifile = 1 : length(fnames)
    if ~fnames(ifile).isdir
        filename = strrep(fnames(ifile).name, basename, recname);
        newfile = fullfile(recpath, filename);
        movefile(fnames(ifile).name, newfile, 'f')
    end
end

% emg data
mapch = [1, 2];
rmvch = [];
store = 'EMG1';
datInfo = tdt2dat('basepath', basepath, 'store', store, 'blocks',  blocks,...
    'chunksize', chunksize, 'mapch', mapch, 'rmvch', rmvch, 'clip', clip);

% create lfp file
LFPfromDat('basepath', basepath, 'cf', 450, 'chunksize', 5e6,...
    'nchans', 2, 'fsOut', 1250, 'fsIn', datInfo.fs)  

% change filenames
movefile([basename, '.lfp'], [basename, '.emg.dat'], 'f')
delete([basename, '.dat'])

% move files to session folder
fnames = dir(['*' basename '*']);
for ifile = 1 : length(fnames)
    if ~fnames(ifile).isdir
        filename = strrep(fnames(ifile).name, basename, recname);
        newfile = fullfile(recpath, filename);
        movefile(fnames(ifile).name, newfile, 'f')
    end
end

% cp xml
mousepath = fileparts(basepath);
xmlmouse = fullfile(mousepath, [mname, '.xml']);
xmlrec = fullfile(recpath, [recname, '.xml']);
if exist(xmlmouse, 'file')
    copyfile(xmlmouse, xmlrec, 'f')
end

% move to session folder
cd(recpath)
session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'forceDef', true, 'forceL', true, 'saveVar', true);      
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;
basepath = session.general.basePath;
[~, basename] = fileparts(basepath);

% sleep signals
emg = bz_LoadBinary([basename, '.emg.dat'], 'frequency', 1250, 'nChannels', 2,...
    'channels', 1);
sSig = as_prepSig([basename, '.lfp'], double(emg),...
    'eegCh', [13 : 16], 'emgCh', [], 'saveVar', true, 'emgNchans', [],...
    'eegNchans', nchans, 'inspectSig', true, 'forceLoad', true,...
    'eegFs', 1250, 'emgFs', 1250, 'eegCf', [],...
    'emgCf', [10 450], 'fs', 1250);
clear sSig

% calc spec
lfp = mean(double(bz_LoadBinary([basename, '.lfp'], 'frequency', 1250, 'nChannels', 16,...
    'channels', [13 : 16])), 2);
eeg = double(bz_LoadBinary([basename, '.emg.dat'], 'frequency', 1250, 'nChannels', 2,...
    'channels', 2));
tlfp = [1 : length(lfp)] / 1250;
eeg = [interp1([1 : length(eeg)] / 1250, eeg, tlfp,...
    'spline')]';
spec = calc_spec('sig', [lfp, eeg], 'fs', 1250, 'graphics', true, 'saveVar', true,...
    'padfft', -1, 'winstep', 5, 'logfreq', true, 'ftarget', [],...
    'ch', [{1}, {2}], 'force', true);

% spike detection from temp_wh
[spktimes, ~] = spktimesWh('basepath', basepath, 'fs', fs, 'nchans', nchans,...
    'spkgrp', spkgrp, 'saveVar', true, 'saveWh', true,...
    'graphics', false, 'force', true, 'winWh', [0 Inf]);

% spike rate per tetrode. note that using firingRate requires
% special care becasue spktimes is given in samples and not seconds
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

delete('*temp_wh*')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lh107
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% user input
mname = 'lh107';
basepath = 'I:\Data\lh106\lh106_220525_091000';
blocks = [1 : 2];

% tank to dat (tetrodes)
store = 'Raw2';
mapch = 1 : 16;
rmvch = [];
chunksize = 300;
clip = cell(1, 1);
datInfo = tdt2dat('basepath', basepath, 'store', store, 'blocks',  blocks,...
    'chunksize', chunksize, 'mapch', mapch, 'rmvch', rmvch, 'clip', clip);

% move dat to session folder
[~, basename] = fileparts(basepath);
recname = strrep(basename, 'lh106', mname);
recpath = fullfile(basepath, recname);
mkdir(recpath)
fnames = dir(['*' basename '*']);
for ifile = 1 : length(fnames)
    if ~fnames(ifile).isdir
        filename = strrep(fnames(ifile).name, basename, recname);
        newfile = fullfile(recpath, filename);
        movefile(fnames(ifile).name, newfile, 'f')
    end
end

% emg data
mapch = 1;
rmvch = [];
store = 'EMG2';
[datInfo] = tdt2dat('basepath', basepath, 'store', store, 'blocks',  blocks,...
    'chunksize', chunksize, 'mapch', mapch, 'rmvch', rmvch, 'clip', clip);

% create lfp file
LFPfromDat('basepath', basepath, 'cf', 450, 'chunksize', 5e6,...
    'nchans', 1, 'fsOut', 1250, 'fsIn', datInfo.fs)  

% change filenames
movefile([basename, '.lfp'], [basename, '.emg.dat'], 'f')
delete([basename, '.dat'])

% move files to session folder
fnames = dir(['*' basename '*']);
for ifile = 1 : length(fnames)
    if ~fnames(ifile).isdir
        filename = strrep(fnames(ifile).name, basename, recname);
        newfile = fullfile(recpath, filename);
        movefile(fnames(ifile).name, newfile, 'f')
    end
end

% cp xml
mousepath = fileparts(basepath);
xmlmouse = fullfile(mousepath, [mname, '.xml']);
xmlrec = fullfile(recpath, [recname, '.xml']);
if exist(xmlmouse, 'file')
    copyfile(xmlmouse, xmlrec, 'f')
end

% move to session folder
cd(recpath)
session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'forceDef', true, 'forceL', true, 'saveVar', true);      
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;
basepath = session.general.basePath;
[~, basename] = fileparts(basepath);

% sleep signals
emg = bz_LoadBinary([basename, '.emg.dat'], 'frequency', 1250, 'nChannels', 1,...
    'channels', 1);
sSig = as_prepSig([basename, '.lfp'], double(emg),...
    'eegCh', [5 : 8], 'emgCh', [], 'saveVar', true, 'emgNchans', [],...
    'eegNchans', nchans, 'inspectSig', true, 'forceLoad', true,...
    'eegFs', 1250, 'emgFs', 1250, 'eegCf', [],...
    'emgCf', [10 450], 'fs', 1250);
clear sSig

% calc spec
spec = calc_spec('sig', [], 'fs', 1250, 'graphics', true, 'saveVar', true,...
    'padfft', -1, 'winstep', 5, 'logfreq', true, 'ftarget', [],...
    'ch', [{1 : 4}], 'force', true);

% spike detection from temp_wh
[spktimes, ~] = spktimesWh('basepath', basepath, 'fs', fs, 'nchans', nchans,...
    'spkgrp', spkgrp, 'saveVar', true, 'saveWh', true,...
    'graphics', false, 'force', true, 'winWh', [0 Inf]);

% spike rate per tetrode. note that using firingRate requires
% special care becasue spktimes is given in samples and not seconds
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

delete('*temp_wh*')




