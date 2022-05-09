    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lh101
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% user input
mname = 'lh101';
basepath = 'G:\Data\lh104\lh104_220503_091700';
blocks = [1, 10, 11, 20, 21];

% tank to dat (tetrodes)
store = 'Raw3';
mapch = 1 : 16;
rmvch = [];
chunksize = 300;
clip = cell(1, 1);
datInfo = tdt2dat('basepath', basepath, 'store', store, 'blocks',  blocks,...
    'chunksize', chunksize, 'mapch', mapch, 'rmvch', rmvch, 'clip', clip);

% move dat to session folder
[~, basename] = fileparts(basepath);
recname = strrep(basename, 'lh104', mname);
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
mapch = 1 : 8;
rmvch = [1 : 7];
store = 'Raw1';
[datInfo, emg] = tdt2dat('basepath', basepath, 'store', store, 'blocks',  blocks,...
    'chunksize', chunksize, 'mapch', mapch, 'rmvch', rmvch, 'clip', clip);

% eeg data
store = 'Raw2';
[datInfo] = tdt2dat('basepath', basepath, 'store', store, 'blocks',  blocks,...
    'chunksize', chunksize, 'mapch', mapch, 'rmvch', rmvch, 'clip', clip);

% change filenames
delete(['*', basename, '.dat'])
movefile([basename, '.Raw1.data.mat'], [basename, '.Raw1.emg.mat'], 'f')
movefile([basename, '.lfp'], [basename, '.eeg'], 'f')

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

% calc spec
spec = calc_spec('sig', [], 'fs', 1250, 'graphics', true, 'saveVar', true,...
    'padfft', -1, 'winstep', 5, 'logfreq', true, 'ftarget', [],...
    'ch', [{13 : 16}], 'force', true);

% sleep signals
clear spec
sSig = as_prepSig([basename, '.lfp'], double(emg),...
    'eegCh', [13 : 16], 'emgCh', [], 'saveVar', true, 'emgNchans', [],...
    'eegNchans', nchans, 'inspectSig', true, 'forceLoad', true,...
    'eegFs', 1250, 'emgFs', 6103.515625, 'eegCf', [],...
    'emgCf', [10 450], 'fs', 1250);
clear sSig

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
% lh104 (Raw1) or lh105 (Raw2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% user input
mname = 'lh105';
basepath = 'G:\Data\lh104\lh104_220503_091700';
blocks = [1, 10, 11, 20, 21];

% tank to dat (tetrodes)
store = 'Raw2';
mapch = [1 : 8];
rmvch = [1 : 2 : 7, 8];
chunksize = 300;
clip = cell(1, 1);
datInfo = tdt2dat('basepath', basepath, 'store', store, 'blocks',  blocks,...
    'chunksize', chunksize, 'mapch', mapch, 'rmvch', rmvch, 'clip', clip);

% move dat to session folder
[~, basename] = fileparts(basepath);
recname = strrep(basename, 'lh104', mname);
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

% cp xml
mousepath = fileparts(basepath);
xmlmouse = fullfile(mousepath, [mname, '.xml']);
xmlrec = fullfile(recpath, [recname, '.xml']);
if exist(xmlmouse, 'file')
    copyfile(xmlmouse, xmlrec, 'f')
end

% move to session
cd(recpath)
session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'forceDef', true, 'forceL', true, 'saveVar', true);      
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
basepath = session.general.basePath;
[~, basename] = fileparts(basepath);

% calc spec
spec = calc_spec('sig', [], 'fs', 1250, 'graphics', true, 'saveVar', true,...
    'padfft', -1, 'winstep', 1, 'logfreq', true, 'ftarget', [],...
    'ch', [{1}, {2}], 'force', true);

% sleep signals
sSig = as_prepSig([basename, '.lfp'], [basename, '.lfp'],...
    'eegCh', [1], 'emgCh', [3], 'saveVar', true, 'emgNchans', [],...
    'eegNchans', nchans, 'inspectSig', true, 'forceLoad', true,...
    'eegFs', 1250, 'emgFs', [], 'eegCf', [],...
    'emgCf', [10 450], 'fs', 1250);


ch = 2;
clear axh
fh = figure;
th = tiledlayout(2, 1, 'TileSpacing', 'Compact');
axh(1) = nexttile;
plot_spec(spec, 'ch', ch, 'logfreq', true, 'saveFig', true,...
    'axh', axh)

axh(2) = nexttile;
yval = movmean(squeeze(spec.bands.db(:, :, ch)), 100,  2);
plot(spec.tstamps / 60 / 60,...
    yval(2 : end, :) ./ yval(1, :), 'LineWidth', 2)
legend(spec.bands.bandNames(2 : end))
ylabel('Norm. Spectral Power [dB]')
xlabel('Time [hr]')
linkaxes(axh, 'x')

tidx = cumsum(datInfo.nsec) / 60 / 60;
hold on
plot([tidx; tidx], ylim, '--k')


% mean trace per intensity
    prism = cellfun(@(x) mean(x, 2, 'omitnan'), traces, 'uni', false);
    prism = cell2mat(prism);