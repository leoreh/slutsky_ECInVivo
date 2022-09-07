% wrapper for batch processing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = 'lh111';
forceL = true;
forceA = true;

pcond = ["tempflag"];
ncond = [""];

% load vars from each session
varsFile = ["fr"; "sr"; "spikes"; "st_metrics"; "swv_metrics";...
    "cell_metrics"; "sleep_states"; "ripp.mat"; "datInfo"; "session";...
    "psd_bins"];
varsName = ["fr"; "sr"; "spikes"; "st"; "swv"; "cm"; "ss"; "ripp";...
    "datInfo"; "session"; "psdBins"];
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';

[v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', pcond, 'ncond', ncond,...
    'xlsname', xlsname);
nfiles = length(basepaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% templateCal = ss.info.calibrationData;

for isession = 1 : nfiles

    % file
    basepath = basepaths{isession};
    cd(basepath)
    [~, basename] = fileparts(basepath);

    % print progress
    fprintf('sessions_wrapper: working on session %d of %d, %s\n',...
        isession, nfiles, basename)

    % params
    session = v(isession).session;
    nchans = session.extracellular.nChannels;
    fs = session.extracellular.sr;
    spkgrp = session.extracellular.spikeGroups.channels;


    %     % add timebins to datInfo
    %     [timebins, timepnt] = metaInfo_timebins('reqPnt', [], 'nbins', 2);
    %     winCalc = mat2cell(timebins, [1, 1], 2)';

%     % sleep signals
%     load([basename, '.sleep_sig.mat'], 'info');
%     if ~any(info.eegCh ==7)
%         load([basename, '.acceleration.mat'])
%         sSig = as_prepSig([basename, '.lfp'], acc.mag,...
%             'eegCh', [7 : 10], 'emgCh', [], 'saveVar', true, 'emgNchans', [],...
%             'eegNchans', nchans, 'inspectSig', true, 'forceLoad', true,...
%             'eegFs', 1250, 'emgFs', 1250, 'eegCf', [], 'emgCf', [10 450], 'fs', 1250);
%     else
%         sSig = load([basename, '.sleep_sig.mat']);
%     end
% 
%     % classify
%     ss = as_classify(sSig, 'basepath', basepath, 'inspectLabels', false,...
%         'saveVar', true, 'forceA', true, 'netfile', [],...
%         'graphics', true, 'calData', templateCal);

        % calc psd in states
        psd = psd_states('basepath', basepath, 'sstates', [1, 4],...
            'ch', [7], 'fs', [], 'wins', [0, Inf], 'saveVar', true,...
            'graphics', true, 'forceA', true);

    %
    %     % calc spec
    %     spec = calc_spec('sig', [], 'fs', 1250, 'graphics', true,...
    %         'saveVar', true, 'padfft', -1, 'winstep', 5,...
    %         'ftarget', [], 'ch', {[7 : 10]},...
    %         'force', true, 'logfreq', true);

end

cell_metrics = CellExplorer('basepaths', basepaths);


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
sstates = psd.info.input.sstates;
faxis = v(1).psd.info.input.faxis;
lim_fAxis = faxis >= 0;

fh = figure;
th = tiledlayout(1, length(sstates), 'TileSpacing', 'Compact');
alphaIdx = linspace(0.5, 1, nfiles);

for istate = 1 : nstates
    axh = nexttile;

    ydata = squeeze(psd.psd(istate, :, :));
%     ydata = ydata ./ sum(ydata(lim_fAxis, :), 1);

    ph = plot(faxis, ydata, 'LineWidth', 2);
    for ifile = 1 : nfiles
        ph(ifile).Color(istate) = cfg.colors{sstates(istate)}(istate) - ifile * 0.01;
        ph(ifile).Color(4) = alphaIdx(ifile);
    end
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
clear psd
for ifile = 1 : nfiles
    
    cd(basepaths{ifile})
    spec(ifile) = calc_spec('basepath', basepaths{ifile},...
        'sig', [], 'fs', 1250, 'graphics', false,...
        'saveVar', true, 'padfft', -1, 'winstep', 1,...
        'ftarget', [], 'ch', {[7]},...
        'force', true, 'logfreq', true);
    
    for istate = 1 : length(sstates)
        psd(ifile, istate, :) = mean(spec(ifile).s(v(ifile).ss.labels == sstates(istate), :), 1);
    end
end

fh = figure;
th = tiledlayout(1, length(sstates), 'TileSpacing', 'Compact');
alphaIdx = linspace(0.5, 1, nfiles);

for istate = 1 : nstates
    axh = nexttile;

    ydata = squeeze(psd(:, istate, :));
%     ydata = ydata ./ sum(ydata(lim_fAxis, :), 1);

    ph = plot(spec(ifile).freq, ydata, 'LineWidth', 2);
    for ifile = 1 : nfiles
        ph(ifile).Color(istate) = cfg.colors{sstates(istate)}(istate) - ifile * 0.01;
        ph(ifile).Color(4) = alphaIdx(ifile);
    end
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
for isession = 1 : nfiles
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



