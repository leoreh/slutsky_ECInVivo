% wrapper for batch processing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = 'lh96';
forceL = true;
forceA = true;

pcond = ["tempflag"];
ncond = [""];

% load vars from each session
varsFile = ["fr"; "sr"; "spikes"; "st_metrics"; "swv_metrics";...
    "cell_metrics"; "sleep_states"; "ripp.mat"; "datInfo"; "session"];
varsName = ["fr"; "sr"; "spikes"; "st"; "swv"; "cm"; "ss"; "ripp";...
    "datInfo"; "session"];
if ~exist('v', 'var') || forceL
    [v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', pcond, 'ncond', ncond);
end
nsessions = length(basepaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
templateCal = ss.info.calibrationData;

if forceA
    for isession = 2 : nsessions
        
        % file
        basepath = basepaths{isession};
        cd(basepath)
        [~, basename] = fileparts(basepath);
        
        % params
        session = v(isession).session;
        nchans = session.extracellular.nChannels;
        fs = session.extracellular.sr;
        spkgrp = session.extracellular.spikeGroups.channels;
        
        % update units
        %         units = selectUnits('basepath', pwd, 'grp', [], 'saveVar', true,...
        %             'forceA', true, 'frBoundries', [0.1 Inf; 0.1 Inf],...
        %             'spikes', []);
        
        % get ripples
%         ripp = getRipples('basepath', basepath, 'rippCh', [9],...
%             'emg', acc.mag, 'recWin', [0, Inf], 'saveVar', true,...
%             'spkFlag', true, 'graphics', true, 'saveVar', true);
        
        %         swv = spkwvMetrics('basepath', basepath, 'fs', fs, 'forceA', true);
        
        
        
        %         % create emg signal from accelerometer data
        %         acc = EMGfromACC('basepath', basepath, 'fname', [basename, '.lfp'],...
        %             'nchans', nchans, 'ch', nchans - 2 : nchans, 'saveVar', true, 'fsIn', 1250,...
        %             'graphics', false, 'force', true);
        %
        %         % state
        load([basename, '.acceleration.mat'])
                sSig = as_prepSig([basename, '.lfp'], acc.mag,...
                    'eegCh', [16], 'emgCh', [], 'saveVar', true, 'emgNchans', [],...
                    'eegNchans', nchans, 'inspectSig', false, 'forceLoad', true,...
                    'eegFs', 1250, 'emgFs', 1250, 'eegCf', [], 'emgCf', [10 450], 'fs', 1250);
        %
        %         % classify with a network
%                 calData = templateCal;
%                 netfile = [];
%                 ss = as_classify(sSig, 'basepath', pwd, 'inspectLabels', false,...
%                     'saveVar', true, 'forceA', true, 'netfile', netfile,...
%                     'graphics', true, 'calData', calData);
        %
        
        %         sSig = load([basename, '.sleep_sig.mat']);
        %         load([basename, '.sleep_states.mat']);
        %         as_stateSeparation(sSig, ss)
        
        
        %         tbins_txt = {'0-3ZT', '3-6ZT', '6-9ZT', '9-12ZT',...
        %             '12-15ZT', '15-18ZT', '18-21ZT', '21-24ZT'};
        %         psdBins = psd_states_timebins('basepath', pwd,...
        %             'chEeg', [], 'forceA', true, 'graphics', true,...
        %             'timebins', chunks, 'saveVar', true,...
        %             'sstates', [1, 4, 5], 'tbins_txt', tbins_txt);
        
        %         psdBins = psd_timebins('basepath', pwd,...
        %             'chEeg', [34], 'forceA', true, 'graphics', true,...
        %             'timebins', chunks, 'saveVar', true);
        
        %         fr = firingRate(v(isession).spikes.times, 'basepath', basepath, 'graphics', true,...
        %             'binsize', 60, 'saveVar', true, 'smet', 'GK', 'winBL',...
        %             [0 timepoints], 'winCalc', [0, Inf], 'forceA', true);
        
        %         frBins(isession) = fr_timebins('basepath', pwd,...
        %             'forceA', false, 'graphics', true,...
        %             'timebins', chunks, 'saveVar', true);
        
    end
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
% states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% rearrange state of entire experiment
sessionIdx = 1 : nsessions;
stateidx = [1, 4, 5];
ts = 1;                         % state labels epoch length
lightHr = '080000';             % when does light cycle start [HHMM]

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

