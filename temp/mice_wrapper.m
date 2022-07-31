   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lh110
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% user input
mname           = 'lh110';
basepath        = 'G:\Data\lh110\lh110_220731_102200';
blocks          = [1 : 2];

% tank to dat 
mapch           = [];
rmvch           = [1, 3];
store           = 'Raw2';
chunksize       = 300;
clip            = cell(1, 1);
datInfo = tdt2dat('basepath', basepath, 'store', store, 'blocks',  blocks,...
    'chunksize', chunksize, 'mapch', mapch, 'rmvch', rmvch, 'clip', clip);

% fepsp 
intens          = [300, 500, 700];
protocol_id     = 'pair';
ch              = 2;
blocks          = 2;                            
fepsp_tdtPipeline('basepath', basepath, 'blocks', blocks,...
    'protocol_id', protocol_id, 'recsuffix', '', 'intens', intens,...
    'ch', ch, 'mapch', mapch', 'rmvch', rmvch, 'store', store)

% move files to session folder
[~, basename] = fileparts(basepath);
recpath = fullfile(basepath, basename);
mkdir(recpath)
fnames = dir(['*' basename '*']);
for ifile = 1 : length(fnames)
    if ~fnames(ifile).isdir
        filename = fnames(ifile).name;
        newfile = fullfile(recpath, filename);
        movefile(fnames(ifile).name, newfile, 'f')
    end
end
movefile('graphics', recpath, 'f')

% move xml file
mousepath = fileparts(basepath);
xmlfile = dir(fullfile(mousepath, [mname, '.xml']));
newname = strrep(xmlfile.name, mname, basename);
newfile = fullfile(recpath, newname);
copyfile(fullfile(mousepath, xmlfile.name), newfile)

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
sSig = as_prepSig([basename, '.lfp'], [],...
    'eegCh', [1], 'emgCh', [2], 'saveVar', true, 'emgNchans', nchans,...
    'eegNchans', nchans, 'inspectSig', true, 'forceLoad', true,...
    'eegFs', 1250, 'emgFs', 1250, 'eegCf', [], 'emgCf', [50 450], 'fs', 1250);

% calc spec
spec = calc_spec('sig', [], 'fs', 1250, 'graphics', true, 'saveVar', true,...
    'padfft', -1, 'winstep', 5, 'logfreq', true, 'ftarget', [],...
    'ch', [{1}], 'force', true);

% -------------------------------------------------------------------------
% fepsp in relation to states

% load data
varsFile = ["fepsp_traces"; "fepsp_results"; "sleep_states"];
varsName = ["traces"; "results"; "ss"];
v = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
    'varsName', varsName);

% get stim indices in specific state
istate = 1;
stateStim = InIntervals(v.results.info.stimIdx, v.ss.stateEpochs{istate});

nstims = 45;
nintens = length(intens);
cellIdx = num2cell(reshape(1 : nstims, nstims / nintens, nintens), 1);
clear amp traces
for iintens = 1 : nintens
    stimIdx = stateStim(cellIdx{iintens});
    amp{iintens} = v.results.all_traces.Amp{iintens}(:, stimIdx);
    traces(iintens, :) = mean(v.traces{iintens}(:, stimIdx), 2);
end

protocol_info = fepsp_getProtocol("protocol_id", 'pair', "fs", fs);
protocol_info.Tstamps;

cell2nanmat(amp, 2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lh109
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% user input
mname           = 'lh109';
basepath        = 'F:\Data\lh110\lh110_220729_111000';
blocks          = [1];
cd(basepath)

% tank to dat 
mapch           = [];
rmvch           = [1, 3, 5];
store           = 'Raw1';
chunksize       = 300;
clip            = cell(1, 1);
datInfo = tdt2dat('basepath', basepath, 'store', store, 'blocks',  blocks,...
    'chunksize', chunksize, 'mapch', mapch, 'rmvch', rmvch, 'clip', clip);

% move files to session folder
[~, basename] = fileparts(basepath);
recname = strrep(basename, 'lh110', mname);
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

% move xml file
mousepath = fileparts(basepath);
xmlfile = dir(fullfile(mousepath, [mname, '.xml']));
newname = strrep(xmlfile.name, mname, recname);
newfile = fullfile(recpath, newname);
copyfile(fullfile(mousepath, xmlfile.name), newfile)

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
sSig = as_prepSig([basename, '.lfp'], [],...
    'eegCh', [1], 'emgCh', [2], 'saveVar', true, 'emgNchans', nchans,...
    'eegNchans', nchans, 'inspectSig', true, 'forceLoad', true,...
    'eegFs', 1250, 'emgFs', 1250, 'eegCf', [], 'emgCf', [50 450], 'fs', 1250);

% calc spec
spec = calc_spec('sig', [], 'fs', 1250, 'graphics', true, 'saveVar', true,...
    'padfft', -1, 'winstep', 5, 'logfreq', true, 'ftarget', [],...
    'ch', [{1}, {3}], 'force', true);

% plot
plot_spec(spec, 'ch', 2, 'logfreq', true, 'saveFig', false,...
    'axh', [])


