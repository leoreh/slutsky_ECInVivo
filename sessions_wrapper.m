% wrapper for batch processing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = 'lh87';
forceL = true;
forceA = true;

pcond = ["tempflag"];
ncond = [""];

% load vars from each session
varsFile = ["fr"; "sr"; "spikes"; "st_metrics"; "swv_metrics";...
    "cell_metrics"; "sleep_states"; "ripp.mat"; "datInfo"; "session"];
varsName = ["fr"; "sr"; "spikes"; "st"; "swv"; "cm"; "ss"; "ripp";...
    "datInfo"; "session"];
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';

if ~exist('v', 'var') || forceL
    [v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', pcond, 'ncond', ncond,...
        'xlsname', xlsname);
end
nsessions = length(basepaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
templateCal = ss.info.calibrationData;

for isession = 17 : nsessions
    
    % file
    basepath = basepaths{isession};
    cd(basepath)
    [~, basename] = fileparts(basepath);

    % print progress
    fprintf('sessions_wrapper: working on session %d of %d, %s\n',...
        isession, nsessions, basename)

    % params
    session = CE_sessionTemplate(pwd, 'viaGUI', false,...
        'forceDef', true, 'forceL', false, 'saveVar', true);
    %     session = v(isession).session;
    nchans = session.extracellular.nChannels;
    fs = session.extracellular.sr;
    spkgrp = session.extracellular.spikeGroups.channels;
    
    % add timebins to datInfo
    [timebins, timepnt] = metaInfo_timebins('reqPnt', [], 'nbins', 2);
    winCalc = mat2cell(timebins, [1, 1], 2)';
    
%     % sleep signals
%     sSig = as_prepSig([basename, '.lfp'], [],...
%         'eegCh', [5], 'emgCh', [3], 'saveVar', true, 'emgNchans', [],...
%         'eegNchans', nchans, 'inspectSig', false, 'forceLoad', true,...
%         'eegFs', 1250, 'emgFs', [], 'eegCf', [], 'emgCf', [10 450], 'fs', 1250);
%     
%     % manual lables
%     labelsfile = [basename, '.AccuSleep_labelsMan.mat'];
%     labelsmanfile = [basename, '.sleep_labelsMan.mat'];
%     if exist(labelsfile, 'file')
%         load(labelsfile, 'labels')
%     end
%     save(labelsmanfile, 'labels')
%     %     AccuSleep_viewer(sSig, labels, labelsmanfile)
%     
%     % classify
%     ss = as_classify(sSig, 'basepath', pwd, 'inspectLabels', false,...
%         'saveVar', true, 'forceA', true, 'netfile', [],...
%         'graphics', true, 'calData', templateCal);

    % calc psd in states
    psdBins = psd_states_timebins('basepath', pwd,...
        'chEeg', [], 'forceA', true, 'graphics', true,...
        'timebins', timebins, 'saveVar', true, 'sstates', [1, 4]);

%     % load data
%     sig = double(bz_LoadBinary([basename, '.lfp'],...
%         'duration', Inf, 'frequency', 1250, 'nchannels', nchans,...
%         'start', 0, 'channels', 2, 'downsample', 1));
%     
%     % calc psd
%     freq = [0.5 : 0.5 : 120];
%     psd = psd_timebins('sig', sig, 'winCalc', winCalc,...
%         'fs', 1250, 'graphics', true, 'faxis', freq);
% 
%     % calc spec
%     spec = calc_spec('sig', [], 'fs', 1250, 'graphics', true,...
%         'saveVar', true, 'padfft', -1, 'winstep', 5,...
%         'ftarget', [], 'ch', [{1}, {2}, {4}, {5}],...
%         'force', true, 'logfreq', true);

end


cell_metrics = CellExplorer('basepaths', basepaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'force', true, 'saveVar', true);
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;

sessionIdx = 1 : nsessions;
stateidx = [1, 4, 5];
grp = [1 : 4];                  % which tetrodes to plot
unitClass = 'pyr';              % plot 'int', 'pyr', or 'all'
suFlag = 1;                     % plot only su or all units
frBoundries = [0 Inf];          % include only units with fr greater than

[nsub] = numSubplots(length(sessionIdx));
[cfg_colors, cfg_names, ~] = as_loadConfig([]);
setMatlabGraphics(false)

% arrange title names
for isession = 1 : nsessions
    sessionName{isession} = dirnames{isession}(length(mname) + 2 : end);
    basepath = char(fullfile(mousepath, dirnames{isession}));
    basepaths{isession} = fullfile(mousepath, dirnames{isession});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% burstiness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load vars from each session
varsFile = ["datInfo"; "session"; "st_metrics"; "units"];
varsName = ["datInfo"; "session"; "st"; "units"];
[v, ~] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);
nsessions = length(basepaths);

% brst vars to organize and plot
% from avg sub-struct
brstVar = {'freq', 'dur', 'spkprct', 'ibi', 'nspks'};

% fron nbrsts sub-struct

unitType = 'rs';

% initialize
brst = cell2struct(cell(1, length(brstVar)), brstVar, 2);

% concate to cell
for isession = 1 : nsessions
    su = v(isession).units.(unitType);

    for ivar = 1 : length(brstVar)
        brst.(brstVar{ivar}){isession} = v(isession).st.brst.avg.(brstVar{ivar})(:, su)';
    end

    brst.nbrsts{isession} = v(isession).st.brst.nbrsts.freqNorm(:, su)';
end

% organize in mat
for ivar = 1 : length(brstVar)
    brst.(brstVar{ivar}) = cell2nanmat(brst.(brstVar{ivar}), 1);
end
brst.nbrsts = cell2nanmat(brst.nbrsts, 1);

% xlabels
xdata = [-5 * 6 : 6 : 136];

% graphics
setMatlabGraphics(true)
fh = figure;
th = tiledlayout(length(brstVar) + 1, 1, 'TileSpacing', 'Compact');

for ivar = 1 : length(brstVar)
    nexttile
    plot_boxMean('dataMat', brst.(brstVar{ivar}), 'clr', 'k', 'allPnts', true)
    xticklabels(xdata);
    xlabel('Time [h]')
    ylabel(brstVar{ivar})
end

nexttile
plot_boxMean('dataMat', brst.nbrsts, 'clr', 'k', 'allPnts', true)
xticklabels(xdata);
xlabel('Time [h]')
ylabel('Norm. Freq')



