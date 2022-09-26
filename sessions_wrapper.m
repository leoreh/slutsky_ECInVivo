% wrapper for batch processing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = 'lh96';

varsFile = ["fr"; "sr"; "spikes"; "st_metrics"; "swv_metrics";...
    "cell_metrics"; "sleep_states"; "ripp.mat"; "datInfo"; "session";...
    "psd_bins"; "units"];
varsName = ["fr"; "sr"; "spikes"; "st"; "swv"; "cm"; "ss"; "ripp";...
    "datInfo"; "session"; "psdBins"; "units"];
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
[v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
    'xlsname', xlsname);
nfiles = length(basepaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
templateCal = ss.info.calibrationData;

for ifile = 1 : nfiles

    % file
    basepath = basepaths{ifile};
    cd(basepath)
    [~, basename] = fileparts(basepath);

    % print progress
    fprintf('sessions_wrapper: working on session %d of %d, %s\n',...
        ifile, nfiles, basename)

    % params
    session = v(ifile).session;
    nchans = session.extracellular.nChannels;
    fs = session.extracellular.sr;
    spkgrp = session.extracellular.spikeGroups.channels;
    
    % spike timing metrics
    st = spktimes_metrics('spikes', v(ifile).spikes, 'sunits', [],...
        'bins', [0 Inf], 'forceA', true, 'saveVar', true, 'fullA', false);

    % select specific units
%     units = selectUnits('basepath', basepath, 'grp', [1 : 4], 'saveVar', true,...
%         'forceA', true, 'frBoundries', [0.05 Inf; 0.05 Inf],...
%         'spikes', v(ifile).spikes);
    
% fr
% fr = firingRate(v(ifile).spikes.times, 'basepath', basepath,...
%     'graphics', true, 'binsize', 60, 'saveVar', true,...
%     'smet', 'none', 'winBL', [0 Inf], 'winCalc', [0, Inf]);

    % sleep signals
%     sSig = as_prepSig([basename, '.lfp'], [],...
%         'eegCh', [2], 'emgCh', [2], 'saveVar', true, 'emgNchans', nchans,...
%         'eegNchans', nchans, 'inspectSig', false, 'forceLoad', true,...
%         'eegFs', 1250, 'emgFs', 1250, 'eegCf', [], 'emgCf', [100 450],...
%         'fs', 1250, 'sigWin', [0 Inf]);

    % classify
%     ss = as_classify(sSig, 'basepath', basepath, 'inspectLabels', false,...
%         'saveVar', true, 'forceA', true, 'netfile', [],...
%         'graphics', true, 'calData', templateCal);

    % calc psd in states
%     psd = psd_states('basepath', basepath, 'sstates', [1, 4],...
%         'ch', [2], 'fs', [], 'wins', [0, Inf], 'saveVar', true,...
%         'graphics', false, 'forceA', true, 'ftarget', [0.1 : 0.5 : 100]);

    % calc spec
%     spec = calc_spec('sig', [], 'fs', 1250, 'graphics', true,...
%         'saveVar', true, 'padfft', -1, 'winstep', 5,...
%         'ftarget', [], 'ch', {[1 : 4], [5 : 8], [13 : 16], [17 : 20], [21 : 24], 25},...
%         'force', true, 'logfreq', true);

end

cell_metrics = CellExplorer('basepaths', basepaths);

% concatenate var from different sessions
[expData, xData] = sessions_catVarTime('mname', mname,...
    'dataPreset', {'fr'}, 'graphics', true,...
    'basepaths', {}, 'xTicksBinsize', 6, 'markRecTrans', true);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psd across sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load vars from each session
varsFile = ["session"; "psd"];
varsName = ["session"; "psd"];
[v, ~] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);
nfiles = length(basepaths);

% cat psd struct
psd = catfields([v(1 : nfiles).psd], 'catdef', 'addim');

% get params
cfg = as_loadConfig();
nstates = size(psd.psd, 1);
sstates = psd.info.sstates(:, 1);
faxis = v(1).psd.info.faxis;
lim_fAxis = faxis >= 0 & faxis < Inf;

fh = figure;
th = tiledlayout(1, length(sstates), 'TileSpacing', 'Compact');
alphaIdx = linspace(0.5, 1, nfiles);

for istate = 1 : nstates
    axh = nexttile;

    ydata = squeeze(psd.psd(istate, :, :));
    ydata = ydata ./ sum(ydata(lim_fAxis, :), 1);
    
%     for ifile = 1 : nfiles
%         ydata(:, ifile) = bz_NormToRange(ydata(:, ifile), [0 1]);
%     end

    ph = plot(faxis, ydata, 'LineWidth', 2);
%     for ifile = 1 : nfiles
%         ph(ifile).Color(istate) = cfg.colors{sstates(istate)}(istate) - ifile * 0.01;
%         ph(ifile).Color(4) = alphaIdx(ifile);
%     end
    set(gca, 'YScale', 'log', 'XScale', 'log')
    title(cfg.names{sstates(istate)})
    xlabel('Frequency [Hz]')
    legend(split(num2str(1 : nfiles)), 'Location', 'Southwest', 'FontSize', 9)

end

% ----------------

% load vars from each session
varsFile = ["session"; "spec"; "sleep_states"];
varsName = ["session"; "spec"; "ss"];
[v, ~] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);
nfiles = length(basepaths);

sstates = [1, 4];
clear psd spec
for ifile = 1 : nfiles

    cd(basepaths{ifile})
    basepath = basepaths{ifile};
    [~, basename] = fileparts(basepath);
    nchans = v(ifile).session.extracellular.nChannels;
    ch = 7;
    dur = 6 * 60 * 60;

    sig = mean(double(bz_LoadBinary([basename, '.lfp'], 'duration', dur,...
        'frequency', 1250, 'nchannels', nchans, 'start', 0,...
        'channels', 7, 'downsample', 1)), 2);

    spec(ifile) = calc_spec('basepath', basepaths{ifile},...
        'sig', sig, 'fs', 1250, 'graphics', false,...
        'saveVar', false, 'padfft', -1, 'winstep', 1,...
        'ftarget', [0.1 : 0.5 : 100], 'ch', {[1]},...
        'force', true, 'logfreq', true);
    
    for istate = 1 : length(sstates)
        stateIdx = v(ifile).ss.labels == sstates(istate);
        stateIdx = stateIdx(1 : dur);
        psd(ifile, istate, :) = mean(spec(ifile).s(stateIdx, :), 1);
    end
end

fh = figure;
th = tiledlayout(1, length(sstates), 'TileSpacing', 'Compact');
alphaIdx = linspace(0.5, 1, nfiles);

for istate = 1 : nstates
    axh = nexttile;

    ydata = squeeze(psd(:, istate, :))';
    ydata = ydata ./ sum(ydata(lim_fAxis, :), 1);

    ph = plot(spec(ifile).freq, ydata, 'LineWidth', 2);
%     for ifile = 1 : nfiles
%         ph(ifile).Color(istate) = cfg.colors{sstates(istate)}(istate) - ifile * 0.01;
%         ph(ifile).Color(4) = alphaIdx(ifile);
%     end
    set(gca, 'YScale', 'log', 'XScale', 'log')
    title(cfg.names{sstates(istate)})
    xlabel('Frequency [Hz]')
    legend(split(num2str(1 : nfiles)), 'Location', 'Southwest', 'FontSize', 9)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% burstiness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load vars from each session
varsFile = ["datInfo"; "session"; "st_metrics"; "units"];
varsName = ["datInfo"; "session"; "st"; "units"];
[v, ~] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);
nfiles = length(basepaths);

% brst vars to organize and plot
% from avg sub-struct
brstVar = {'freq', 'dur', 'spkprct', 'ibi', 'nspks'};

% fron nbrsts sub-struct

unitType = 'rs';

% initialize
brst = cell2struct(cell(1, length(brstVar)), brstVar, 2);

% concate to cell
for ifile = 1 : nfiles
    su = v(ifile).units.(unitType);

    for ivar = 1 : length(brstVar)
        brst.(brstVar{ivar}){ifile} = v(ifile).st.brst.avg.(brstVar{ivar})(:, su)';
    end

    brst.nbrsts{ifile} = v(ifile).st.brst.nbrsts.freqNorm(:, su)';
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



