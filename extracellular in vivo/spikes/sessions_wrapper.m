% wrapper for batch processing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = 'lh93';
forceL = false;
forceA = false;

% full path and name to xls file with session metadata
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';

% conditions
pcond = ["tempflag"];
ncond = ["fepsp"];

% string array of variables to load
vars = ["session.mat";...
    "cell_metrics.cellinfo";...
    "spikes.cellinfo";...
    "fr.mat";...
    "datInfo";...
    "AccuSleep_states";...
    "sr.mat"];

if ~exist('varArray', 'var') && ~forceL
    [varArray, dirnames, mousepath] = getSessionVars('vars', vars,...
        'pcond', pcond, 'ncond', ncond, 'sortDir', false, 'dirnames', [],...
        'xlsname', xlsname, 'mname', mname);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nsessions = length(dirnames);
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
if length(dirnames) > 1
    pathPieces = regexp(dirnames(:), '_', 'split'); % assumes filename structure: animal_date_time
    sessionDate = [pathPieces{:}];
    sessionDate = sessionDate(2 : 3 : end);
else
    pathPieces = regexp(dirnames(:), '_', 'split');
    sessionDate = {pathPieces{2}};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if forceA
    for isession = 1 : nsessions
        
        % file
        basepath = char(fullfile(mousepath, dirnames{isession}));
        basepaths{isession} = fullfile(mousepath, dirnames{isession});
        cd(basepath)
        [~, basename] = fileparts(basepath);
        
        % params
        session = CE_sessionTemplate(pwd, 'viaGUI', false,...
            'force', true, 'saveVar', true);
        nchans = session.extracellular.nChannels;
        fs = session.extracellular.sr;
        spkgrp = session.extracellular.spikeGroups.channels;
        
        % vars
        % assignVars(varArray, isession)
        
        cell_metrics = ProcessCellMetrics('session', session,...
            'manualAdjustMonoSyn', false, 'summaryFigures', false,...
            'debugMode', true, 'transferFilesFromClusterpath', false,...
            'submitToDatabase', false, 'getWaveformsFromDat', true);
%                 cell_metrics = CellExplorer('basepath', basepath);
        
        load([basename '.spikes.cellinfo.mat'])
        cc = cellclass('basepath', basepath,...
            'waves', cat(1, spikes.rawWaveform{:})', 'saveVar', true,...
            'graphics', false, 'fs', fs);
        
        spikes = cluVal('spikes', spikes, 'basepath', basepath, 'saveVar', true,...
            'saveFig', false, 'force', true, 'mu', [], 'graphics', false,...
            'vis', 'on', 'spkgrp', spkgrp);
        
        %         load([basename, '.AccuSleep_EEG.mat'])
        %         load([basename, '.AccuSleep_EMG.mat'])
        %         %         labelsmanfile = [basename, '.AccuSleep_labelsMan.mat'];
        %         %         AccuSleep_viewer(EEG, EMG, 1250, 1, labels, labelsmanfile)
        %         netfile = 'D:\Code\slutskycode\extracellular in vivo\lfp\SleepStates\AccuSleep\trainedNetworks\net_210708_200155.mat';
        %         ss = as_wrapper(EEG, EMG, [], 'basepath', basepath, 'calfile', [],...
        %             'viaGui', false, 'forceCalibrate', false, 'inspectLabels', false,...
        %             'saveVar', false, 'forceAnalyze', true, 'fs', 1250, 'netfile', [],...
        %             'graphics', true);
        
        % firing rate
%         load([basename '.spikes.cellinfo.mat'])
        binsize = 60;
        winBL = [1 * 60 110 * 60];
        fr = firingRate(spikes.times, 'basepath', basepath,...
            'graphics', false, 'saveFig', false,...
            'binsize', binsize, 'saveVar', true, 'smet', 'MA',...
            'winBL', winBL);
%         
%         %         spike rate per tetrode. note that using firingRate requires
%         %         special care becasue spktimes is given in samples and not seconds
%         load(fullfile(basepath, [basename '.spktimes.mat']))
%         for ii = 1 : length(spkgrp)
%             spktimes{ii} = spktimes{ii} / fs;
%         end
%         binsize = 60;
%         sr = firingRate(spktimes, 'basepath', basepath,...
%             'graphics', false, 'saveFig', false,...
%             'binsize', binsize, 'saveVar', 'sr', 'smet', 'none',...
%             'winBL', [0 Inf]);
    end
end

cell_metrics = CellExplorer('basepaths', basepaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% number of SU and MU (top - pyr; bottom - int) 
grp = [1 : 4];          % which tetrodes to plot
suFlag = 0;             % plot only su or all units
frBoundries = [0 Inf];  % include only units with fr greater / lower than
figFlag = 1; saveFig = true;

if figFlag
    clear units
    for isession = 1 : nsessions
        assignVars(varArray, isession)
        units(isession, 1) = sum(selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'pyr'));
        units(isession, 2) = sum(selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'int'));
    end
    fh = figure;
    bar(units, 'stacked')
    legend({"RS"; "FS"})
    xticks(1 : nsessions)
    xticklabels(dirnames)
    xtickangle(45)
    title('Number of Units')
    xlabel('Session')
    ylabel('No. Units')
    box off
    
    if saveFig
        figname = fullfile(mousepath, 'graphics', 'UnitsDetected');
        export_fig(figname, '-tif', '-transparent', '-r300')
    end
end

% -------------------------------------------------------------------------
% mean firing rate across time for all sessions concatenated (one figure)
grp = [1 : 4];                  % which tetrodes to plot
suFlag = 0;                     % plot only su or all units
frBoundries = [0 Inf];          % include only units with fr greater / lower than
sessionIdx = 1 : nsessions;     % selection of sessions
sessionIdx = 2;

% find date time of start and end of experiment (round to hour)
assignVars(varArray, sessionIdx(end))
ts = sr.binsize;     
[dtStart, ~] = guessDateTime(dirnames(sessionIdx(1)));
dtStart = dateshift(dtStart, 'start', 'hour');
[dtEnd, ~] = guessDateTime(dirnames(sessionIdx(end)));
dtEnd = dtEnd + seconds(length(sr.strd) * ts);
dtEnd = dateshift(dtEnd, 'end', 'hour');
expLen = seconds(dtEnd - dtStart) / ts;     

% get datetime indices and labels for x axis
dtAxis = dtStart : hours(6) : dtEnd;
tidx = 1;
for itime = 1 : length(dtAxis)
    [~, tidx(itime)] = tstamp2time('dtstr', dtStart, 'tstr', dtAxis(itime), 'fs', 1 / ts);
    if tidx(itime) == 0
        tidx(itime) = 1;
    end
    tlabel{itime} = datestr(datenum(dtAxis(itime)), 'dd/mm_HH:MM');
end

% get tidx of manipulation
dt1 = datetime(2021, 02, 28, 19, 00, 00);
dt2 = datetime(2021, 03, 03, 17, 10, 00);
[~, shadeIdx(1)] = tstamp2time('dtstr', dtStart, 'tstr', dt1, 'fs', 1 / ts);
[~, shadeIdx(2)] = tstamp2time('dtstr', dtStart, 'tstr', dt2, 'fs', 1 / ts);
shadeIdx = [0, 0];  % override

% initialize vars that will carry information from entire experiment
if ~isempty(varArray{isession, 3})
    expRS = nan(expLen, length(spikes.su));
    expFS = nan(expLen, length(spikes.su));
end
expMU = nan(expLen, 4);

% concatenate data from different sessions. finds the index of rec start
% from basename and assumes the session recording is continuous
idx_recStart = zeros(length(sessionIdx), 1);
for isession = sessionIdx
    assignVars(varArray, isession)   
    
    % indices
    idx_recStart(isession) = floor(max([1, seconds(guessDateTime(dirnames(isession)) - dtStart) / ts]));
    idx_session = idx_recStart(isession) : idx_recStart(isession) + length(sr.strd) - 1;
    
    % MU
    expMU(idx_session, :) = sr.strd';

    if ~isempty(varArray{isession, 3})
        % RS
        RSunits{isession} = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'pyr');
        expRS(idx_session, 1 : sum(RSunits{isession})) = fr.strd(RSunits{isession}, :)';
        
        % FS
        FSunits{isession} = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'int');
        expFS(idx_session, 1 : sum(FSunits{isession})) = fr.strd(FSunits{isession}, :)';
    end 
end

% plot
fh = figure;
sb1 = subplot(2, 1, 1);
yLimit = [1 max(range(expMU))];
rectangle('Position', [shadeIdx(1) yLimit(1) diff(shadeIdx) diff(yLimit)],...
    'Curvature',0.2, 'faceColor', [0.95 0.95 0.95])
hold on
plot(expMU)
axis tight
xlim([1 length(expMU)])
set(gca, 'box', 'off');
xticks(tidx)
xticklabels(tlabel)
xtickangle(45)
ylabel('Multi-unit firing rate [Hz]')
lgh = legend(split(num2str(grp)));

sb2 = subplot(2, 1, 2);
plot(mean(expRS', 'omitnan'), 'b', 'LineWidth', 2)
hold on
plot(mean(expFS', 'omitnan'), 'r', 'LineWidth', 2)
plot([idx_recStart; idx_recStart], ylim, '--k');
axis tight
xlabel('Time [h]')
ylabel('Single unit firing rate [Hz]')
set(gca, 'box', 'off')
xticks(tidx)
xticklabels(tlabel)
xtickangle(45)
legend(sprintf('RS ~= %d su', round(mean(sum(cell2nanmat(RSunits(sessionIdx), 2), 'omitnan')))),...
 sprintf('FS ~= %d su', round(mean(sum(cell2nanmat(FSunits(sessionIdx), 2), 'omitnan')))));
linkaxes([sb1, sb2], 'x')
if saveFig
    figpath = fullfile(mousepath, 'graphics');
    mkdir(figpath)
    figname = fullfile(figpath, ['FR_time']);
    export_fig(figname, '-tif', '-transparent', '-r300')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% rearrange state of entire experiment
sessionIdx = 1 : nsessions;
stateidx = [1, 4, 5];
ts = 1;                         % state labels epoch length
lightHr = '080000';             % when does light cycle start [HHMM]

[cfg_colors, cfg_names, ~] = as_loadConfig([]);
assignVars(varArray, sessionIdx(end))

% find date time of start and end of experiment and round it according to
% the light cycle
[dtStart, ~] = guessDateTime(dirnames(sessionIdx(1)));
dtStart = dateshift(dtStart, 'start', 'hour');
dtStart = dtStart + timeofday(guessDateTime(lightHr)) - timeofday(dtStart);
[dtEnd, ~] = guessDateTime(dirnames(sessionIdx(end)));
dtEnd = dtEnd + seconds(length(ss.labels) * ts);
dtEnd = dateshift(dtEnd, 'end', 'hour');
dtEnd = dtEnd + timeofday(guessDateTime(lightHr)) - timeofday(dtEnd);
expLen = seconds(dtEnd - dtStart) / ts;     

% initialize vars that will carry information from entire experiment
expLabels = nan(expLen, 1);
expEmg = nan(expLen, 1);

% concatenate data from different sessions. finds the index of rec start
% from basename and assumes the session recording is continuous
for isession = sessionIdx
    assignVars(varArray, isession)   
    basepath = fullfile(mousepath, dirnames{isession});
    cd(basepath)
    [~, basename] = fileparts(basepath);
    idx_recStart = max([1, seconds(guessDateTime(dirnames(isession)) - dtStart) * ts]);

    % labels
    expLabels(idx_recStart : idx_recStart + length(ss.labels) - 1) = ss.labels;   
    
    % emg
    load([basename, '.AccuSleep_EMG.mat'])
    processedEMG = processEMG(standardizeSR(EMG, 1250, 128), 128, 1);  
    processedEMG = (processedEMG - ss.info.calibrationData(end, 1)) ./...
        ss.info.calibrationData(end, 2);
    processedEMG = (processedEMG + 4.5) ./ 9;    
    processedEMG(processedEMG < 0) = 0;
    processedEMG(processedEMG > 1) = 1;
    expEmg(idx_recStart : idx_recStart + length(processedEMG) - 1) = processedEMG;
end

% percent time of state in timebins
dataNorm = 'm';                     % can be '%' or 'm'
binHr = 3;                          % [h]
binSample = binHr * 60 * 60 / ts;   % [s]
bins = n2chunks('n', expLen, 'chunksize', binSample);
nbins = size(bins, 1);
stateBin = nan(nbins, length(stateidx));
for ibin = 1 : nbins
    % data
    for istate = 1 : length(stateidx)
        stateBin(ibin, istate) =...
            sum(expLabels(bins(ibin, 1) : bins(ibin, 2)) ==...
            stateidx(istate), 'omitnan');
        if strcmp(dataNorm, '%')
            stateBin(ibin, istate) = stateBin(ibin, istate) /...
                sum(~isnan(expLabels(bins(ibin, 1) : bins(ibin, 2))));
        elseif strcmp(dataNorm, 'm')
                 stateBin(ibin, istate) = stateBin(ibin, istate) /...
                60 * ts;
        end
    end
    % time labels
    dt = tstamp2time('dtstr', dtStart, 'tstamp', bins(ibin, 1));
    dt = dateshift(dt, 'start', 'hour');
%     tlabel{ibin} = datestr(datenum(dt), 'dd/mm_HH:MM');   
end
tlabel = datestr(datenum(dtStart : hours(binHr) : dtStart + hours(24)), 'HH:MM');
tlegend = datestr(datenum(dtStart : hours(24) : dtEnd), 'dd/mm');

% create 24 hr cycles 
% saveFig = true;
% fh = figure;
% sessionidx = [1, 2];
% clr = ['rrb'];
% for istate = 1 : length(stateidx)
%     subplot(length(stateidx), 1, istate)
%     stateMat{istate} = reshape(stateBin(:, istate), 24 / binHr, size(stateBin(:, istate), 1) / (24 / binHr));
%     stateMat{istate} = stateMat{istate}(:, sessionidx);
%     ph = plot([1.5 : 24 / binHr + 0.5], stateMat{istate}, 'LineWidth', 2);
%     set(ph, {'Color'}, num2cell(clr)')
%     xticks([1 : 24 / binHr + 1])
%     xticklabels(tlabel)
%     xtickangle(45)
%     title(cfg_names(stateidx(istate)))
%     if istate == 1
%         legend(tlegend, 'Location', 'NorthWest')
%     end
%     ylabel(['State Duration [', dataNorm, ']'])
%     xlabel('Time')
% end
% if saveFig
%     figpath = fullfile(mousepath, 'graphics');
%     mkdir(figpath)
%     figname = fullfile(figpath, ['stateDuration_perDay']);
%     export_fig(figname, '-tif', '-transparent', '-r300')
% end
    
% graphics
fh = figure;
subplot(4, 1, 1)
hold on
for istate = stateidx
    stateLabels = find(expLabels == istate);
    scatter(stateLabels / fs / 60 / 60,...
        expEmg(stateLabels),...
        3, cfg_colors{istate})
end
xlim([1 / fs / 60 / 60, length(expLabels) / fs / 60 / 60])
ylim([min(expEmg) 1])
ylabel('Norm. EMG RMS')
set(gca, 'XTick', [], 'YTick', [])

subplot(4, 1, [2 : 4])
xidx = bins(:, 2) - binSample / binHr;
ph = plot(xidx, stateBin, 'Marker', '.', 'LineStyle', '-', 'MarkerSize', 20);
set(ph, {'MarkerEdgeColor'}, cfg_colors(stateidx),...
    {'MarkerFaceColor'}, cfg_colors(stateidx), {'Color'}, cfg_colors(stateidx))
axis tight
hold on
plot([bins(1 : 2 : end, 1), bins(1 : 2 : end, 1)], ylim, '--k', 'LineWidth', 0.5)
xticks(bins(1 : 2 : end, 1))
xticklabels(tlabel(1 : 2 : end, :))
xtickangle(45)
xlabel('Time')
ylabel(['State Duration [', dataNorm, ']'])

if saveFig
    figpath = fullfile(mousepath, 'graphics');
    mkdir(figpath)
    figname = fullfile(figpath, ['stateDuration_aligned']);
    export_fig(figname, '-tif', '-transparent', '-r300')
end

% -------------------------------------------------------------------------
% stacked bar plot of state duration across sessions
figFlag = 1; saveFig = true;
stateidx = [1 : 6];
sessionIdx = 1 : nsessions;
[nsub] = numSubplots(length(sessionIdx));
if figFlag
    clear stateDur
    for isession = sessionIdx
        assignVars(varArray, isession)
        [cfg_colors, cfg_names, ~] = as_loadConfig([]);
        for istate = stateidx 
            stateDur(isession, istate) = sum(ss.epLen{istate}) /...
                length(ss.labels) * 100;
        end               
    end
    
    fh = figure;
    bh = bar(stateDur, 'stacked', 'FaceColor', 'flat');
    for isession = sessionIdx
        for istate = stateidx          
            bh(istate).CData(isession, :) = cfg_colors{istate};
        end
    end
    xticks([1 : length(sessionIdx)]);
    xticklabels(cellstr(dirnames))
    xtickangle(45)
    ylim([0 100])
    ylabel('Time [%]')
    legend(ss.labelNames)
    title('State Duration')
    
    if saveFig
        figpath = fullfile(mousepath, 'graphics');
        mkdir(figpath)
        figname = fullfile(figpath, ['stateDuration']);
        export_fig(figname, '-tif', '-transparent', '-r300')
    end
end

% -------------------------------------------------------------------------
% distribution of mu bins in each state
figFlag = 1; saveFig = true;
grp = 1 : 4;
stateidx = [1 : 6];
sessionIdx = 1 : nsessions;
[nsub] = numSubplots(length(sessionIdx));
plotStyle = 'box';      % can be 'histogram' or 'box'

if figFlag
    fh = figure;
    for isession = 1 : length(sessionIdx)
        
        assignVars(varArray, sessionIdx(isession))
        subplot(nsub(1), nsub(2), isession)
        hold on
        for istate = stateidx
            srmat = sr.states.fr{istate}(grp, :);
            srStates{istate} = srmat(:);
            
            if strcmp(plotStyle, 'histogram')
                binstates{isession, istate} = mean(sr.states.fr{istate}(grp, :), 1);
                h = histogram(binstates{isession, istate}, 20,...
                    'EdgeAlpha', 0, 'faceColor', cfg_colors{istate}, 'FaceAlpha',...
                    0.3, 'Normalization', 'probability');
            end
        end
        
        if strcmp(plotStyle, 'histogram')
            set(gca, 'xscale', 'log', 'box', 'off', 'TickLength', [0 0])
            xlim([10 250])
            ylabel('Counts')
            xlabel('Spike Rate [Hz]')          
            
        elseif strcmp(plotStyle, 'box')
            srmat = cell2nanmat(srStates(stateidx));
            boxplot(srmat, 'PlotStyle', 'traditional', 'Whisker', 6);
            bh = findobj(gca, 'Tag', 'Box');
            bh = flipud(bh);
            for istate = 1 : length(stateidx)
                patch(get(bh(istate), 'XData'), get(bh(istate), 'YData'),...
                    cfg_colors{stateidx(istate)}, 'FaceAlpha', 0.5)
            end
            xticks([])
            ylabel('MU spike rate [Hz]')
            ylim([0 180])
        end
        
        title(sessionDate{sessionIdx(isession)})
        if isession == 1
            legend(ss.labelNames(stateidx))
        end
    end
    
    if saveFig
        figpath = fullfile(mousepath, 'graphics');
        mkdir(figpath)
        figname = fullfile(figpath, ['mu_states']);
        export_fig(figname, '-tif', '-transparent', '-r300')
    end   
end

% -------------------------------------------------------------------------
% distribution of su firing rate per state
figFlag = 1; saveFig = true;
stateidx = [1 : 6];
sessionIdx = 1 : nsessions;
[nsub] = numSubplots(length(sessionIdx));
grp = [1 : 4];          % which tetrodes to plot
unitClass = 'int';      % plot 'int', 'pyr', or 'all'
suFlag = 0;             % plot only su or all units
frBoundries = [0 Inf];  % include only units with mean fr in these boundries

if figFlag
    fh = figure;
    clear frState
    for isession = 1 : length(sessionIdx)
        assignVars(varArray, sessionIdx(isession))
        units = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, unitClass);
        nunits = sum(units);
        
        for istate = stateidx
            frState{isession}(:, istate) = mean(fr.states.fr{istate}(units, :), 2, 'omitnan');
        end
        subplot(nsub(1), nsub(2), isession)
        hold on
        boxplot(frState{isession}(:, stateidx), 'PlotStyle', 'traditional', 'Whisker', 8);
        bh = findobj(gca, 'Tag', 'Box');
        bh = flipud(bh);
        for ibox = 1 : length(bh)
            patch(get(bh(ibox), 'XData'), get(bh(ibox), 'YData'),...
                cfg_colors{stateidx(ibox)}, 'FaceAlpha', 0.5)
        end
        xticks([])
        ylabel([unitClass, ' firing rate [Hz]'])
        ylim([0 15])
        title(sessionDate{sessionIdx(isession)})
        if isession == 1
            legend(ss.labelNames(stateidx))
        end
    end
    if saveFig
        figpath = fullfile(mousepath, 'graphics');
        mkdir(figpath)
        figname = fullfile(figpath, [unitClass, '_states']);
        export_fig(figname, '-tif', '-transparent', '-r300')
    end
end



