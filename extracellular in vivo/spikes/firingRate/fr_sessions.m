% fr_sessions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
forceL = false;
forceA = false;

% full path and name to xls file with session metadata
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
mname = 'lh87';

% column name in xls sheet where dirnames exist
colName = 'Session';

% string array of variables to load
vars = ["session.mat";...
    "cell_metrics.cellinfo";...
    "spikes.cellinfo";...
    "SleepState.states";....
    "fr.mat";...
    "datInfo";...
    "AccuSleep_states";...
    "sr.mat"];

% column name of logical values for each session. only if true than session
% will be loaded. can be a string array and than all conditions must be
% met.
pcond = ["tempflag"; "mancur"; "states"];

% same but imposes a negative condition
ncond = ["fepsp"];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('varArray', 'var') && ~forceL
    [varArray, dirnames, mousepath] = getSessionVars('vars', vars,...
        'pcond', pcond, 'ncond', ncond, 'sortDir', false, 'dirnames', [],...
        'xlsname', xlsname, 'mname', mname);
end
nsessions = length(dirnames);

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
        
        %         % spikes
        %         fixSpkAndRes('grp', [], 'fs', fs, 'stdFactor', 0, 'nchans', nchans);
        %         spikes = loadSpikes('session', session);
        %         spikes = fixCEspikes('basepath', basepath, 'saveVar', false,...
        %             'force', true);
        
%         cell_metrics = ProcessCellMetrics('session', session,...
%             'manualAdjustMonoSyn', false, 'summaryFigures', false,...
%             'debugMode', true, 'transferFilesFromClusterpath', false,...
%             'submitToDatabase', false, 'getWaveformsFromDat', true);
%         cell_metrics = CellExplorer('basepath', basepath);
        
%         load([basename '.spikes.cellinfo.mat'])
%         cc = cellclass('basepath', basepath,...
%             'waves', cat(1, spikes.rawWaveform{:})', 'saveVar', true,...
%             'graphics', false, 'fs', fs);
        
%         spikes = cluVal('spikes', spikes, 'basepath', basepath, 'saveVar', true,...
%             'saveFig', false, 'force', true, 'mu', [], 'graphics', false,...
%             'vis', 'on', 'spkgrp', spkgrp);
        
% load([basename, '.AccuSleep_EEG.mat'])
% load([basename, '.AccuSleep_EMG.mat'])
% ss = as_wrapper(EEG, EMG, [], 'basepath', basepath, 'calfile', [],...
%     'viaGui', false, 'forceCalibrate', false, 'inspectLabels', false,...
%     'saveVar', true, 'forceAnalyze', true, 'fs', 1250, 'netfile', [],...
%     'graphics', true);

        % firing rate
        load([basename '.spikes.cellinfo.mat'])
        binsize = 60;
        winBL = [1 * 60 110 * 60];
        fr = firingRate(spikes.times, 'basepath', basepath,...
            'graphics', false, 'saveFig', false,...
            'binsize', binsize, 'saveVar', true, 'smet', 'MA',...
            'winBL', winBL);
        
%         spike rate per tetrode. note that using firingRate requires
%         special care becasue spktimes is given in samples and not seconds
        load(fullfile(basepath, [basename '.spktimes.mat']))
        for ii = 1 : length(spkgrp)
            spktimes{ii} = spktimes{ii} / fs;
        end
        binsize = 60;
        sr = firingRate(spktimes, 'basepath', basepath,...
            'graphics', false, 'saveFig', false,...
            'binsize', binsize, 'saveVar', 'sr', 'smet', 'none',...
            'winBL', [0 Inf]);               
    end
end

% cell_metrics = CellExplorer('basepaths', basepaths);

% params
session = varArray{1, 1}.session;
session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'force', true, 'saveVar', true);
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% arrange title names
if length(dirnames) > 1
    pathPieces = regexp(dirnames(:), '_', 'split'); % assumes filename structure: animal_date_time
    sessionDate = [pathPieces{:}];
    sessionDate = sessionDate(2 : 3 : end);
else
    pathPieces = regexp(dirnames(:), '_', 'split');
    sessionDate = {pathPieces{2}};
end

setMatlabGraphics(false)

sessionidx = 1 : nsessions;
[nsub] = numSubplots(length(sessionidx));

% close all
grp = [1 : 4];          % which tetrodes to plot
state = [];             % [] - all; 1 - awake; 2 - NREM
FRdata = 'strd';        % plot absolute fr or normalized
unitClass = 'pyr';      % plot 'int', 'pyr', or 'all'
suFlag = 1;             % plot only su or all units
frBoundries = [0 Inf];  % include only units with fr greater than
Y = [0 10];             % ylim
p1 = 1;                 % firing rate vs. time, one fig per session

plotStyle = 'box';      % for p2. can by 'bar' or 'box'
clr = ['bbbkkkkrrr'];


for isession = 1 : nsessions
    
    % vars
    assignVars(varArray, isession)
    fs = session.extracellular.sr;
    
    units{isession} = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, unitClass);
    
    % states
    if ~isempty(ss)
        if isfield(ss, 'ints')
            states = {ss.ints.WAKEstate, ss.ints.NREMstate, ss.ints.REMstate};
        elseif isfield(ss, 'stateEpochs')
            states = ss.stateEpochs;
        end
    else
        states = [];
    end
    if ~isempty(state) && state > 0
        data = fr.states.fr{state}(units{isession}, :);
        tstamps = fr.states.tstamps{state};
        txt = sprintf('mFR of %s %s in %s', unitClass, fr.states.statenames{state});
    else
        data = fr.(FRdata)(units{isession}, :);
        tstamps = fr.tstamps;
        txt = sprintf('mFR of %s %s', unitClass);
    end
    mfr{isession} = mean(data, 2);
    
    if isfield(datInfo, 'nsamps')
        nsamps = cumsum(datInfo.nsamps);
    else
        nsamps = 120 * 60 * fs;
    end
    
    % firing rate vs. time. 1 fig per session
    if p1
        figure
        plot(tstamps / 60, (data'), 'LineWidth', 1)
        hold on
        medata = median((data), 1, 'omitnan');
        plot(tstamps / 60, medata, 'k', 'LineWidth', 4)
        %         stdshade(data, 0.3, 'k', tstamps / 60)
        for ii = 1 : length(nsamps) - 1
            plot([nsamps(ii) nsamps(ii)] / fs / 60, ylim, '--k',...
                'LineWidth', 2)
        end
        axis tight
        ylim(Y)
        %         set(gca, 'YScale', yscale)
        if ~isempty(states)
            fill([states{2} fliplr(states{2})]' / 60, [Y(1) Y(1) Y(2) Y(2)],...
                'b', 'FaceAlpha', 0.15,  'EdgeAlpha', 0);
        end
        ylabel('Firing Rate [Hz]')
        xlabel('Time [m]')
        if isfield(datInfo, 'spktrim')
            xlim(datInfo.spktrim.edges / fs / 60)
        end
        suptitle(dirnames{isession})
        if saveFig
            if nsessions == 1
                clear sessionDate
                sessionDate{1} = pathPieces{2};
            end
            figname = sprintf('%s_%s_FRvsTime', sessionDate{isession}, unitClass)
            suptitle(figname)
            figname = fullfile(mousepath, figname);
            % print(fh, figname, '-dpdf', '-bestfit', '-painters');
            export_fig(figname, '-tif', '-transparent', '-r300')
        end
        
        % fr histogram
        %         fh = figure;
        %         subplot(2, 3, 1)
        %         histogram(log10(mean(fr.states.fr{1}, 2)), 20)
    end
end

figFlag = 0;
if figFlag
    fh = figure;
    mfrmat = cell2nanmat(mfr);
    switch plotStyle
        case 'bar'
            %             plot([1 : nsessions], mfrmat)
            hold on
            bh = bar(nanmean(mfrmat), 'FaceAlpha', 0.1, 'FaceColor', 'k');
            hold on
            % errorbar(1 : nsessions, nanmean(mfrmat), nanstd(mfrmat), 'k')
            idx = [find(clr == 'k', 1) find(clr == 'r', 1)] - 0.5;
            yLimit = ylim;
            patch([idx idx(2) idx(1)],...
                [yLimit(1) yLimit(1) yLimit(2) yLimit(2)],...
                'b', 'FaceAlpha', 0.1)
        case 'box'
            boxplot(mfrmat, 'PlotStyle', 'traditional');
            %             boxplot(mfrmat, 'PlotStyle', 'traditional',...
            %                 'DataLim', [5 10], 'ExtremeMode', 'compress');
            bh = findobj(gca, 'Tag', 'Box');
            bh = flipud(bh);
            if length(bh) == length(clr)
                clrRep = histc(clr, unique(clr));
                clear alphaIdx
                for isession = 1 : length(clrRep)
                    alphaIdx{isession} = linspace(0.8, 0.3, clrRep(isession));
                end
                alphaIdx = [alphaIdx{:}];
                for isession = 1 : length(bh)
                    patch(get(bh(isession), 'XData'), get(bh(isession), 'YData'),...
                        clr(isession), 'FaceAlpha', alphaIdx(isession))
                end
            end
    end
    ylim(Y)
    ylabel('Spike Rate [Hz]')
    xlabel('Session')
    xticks(1 : nsessions)
    xticklabels(sessionDate)
    xtickangle(45)
    box off
    title(txt)
    
    if saveFig
        figname = sprintf('MSR of Tetrode #%s', num2str(grp));
        figname = fullfile(mousepath, figname);
        % print(fh, figname, '-dpdf', '-bestfit', '-painters');
        export_fig(figname, '-tif', '-transparent', '-r300')
    end
end

% number of SU and MU (top - pyr; bottom - int) from selected tetrodes
figFlag = 0;
if figFlag
    clear units
    for isession = 1 : nsessions
        assignVars(varArray, isession)
        units(isession, 1) = sum(selectUnits(spikes, cm, fr, 1, grp, frBoundries, 'pyr'));
        units(isession, 2) = sum(selectUnits(spikes, cm, fr, 1, grp, frBoundries, 'int'));
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

% lfp fft 
figFlag = 0;
ch = [9 : 12];
stateidx = [1, 4];
win = hann(2 ^ nextpow2(10 * fs));
noverlap = 5 * fs;
clear psdPow
if figFlag
    % arrange data
    for isession = 1 : nsessions
        
        assignVars(varArray, isession)
        basepath = session.general.basePath;
        nchans = session.extracellular.nChannels;
        cd(basepath)
        [~, basename] = fileparts(basepath);
        
        % load lfp
        fs = 1250;
        istart = 0;
        dur = 5.9 * 60 * 60;
        lfp = double(bz_LoadBinary([basename, '.lfp'], 'duration', dur,...
            'frequency', 1250, 'nchannels', nchans, 'start', istart,...
            'channels', ch, 'downsample', 1));       
        
        % limit lfp to state
        counter = 1;
        for istate = 1 : length(stateidx)
            epochs = ss.stateEpochs{istate};
            epochInTime = epochs(:, 1) > 2 * 60 * 60 & epochs(:, 2) < dur;
            epochs = epochs(epochInTime, :);
            epochidx = [];
            for iepoch = 1 : size(epochs, 1)
                epochidx = [epochidx, epochs(iepoch, 1) * fs : epochs(iepoch, 2) * fs];
            end
            lfpInState = mean(lfp(epochidx, :), 2);
            lfpInState = lfpInState - mean(lfpInState);
            
            % calc psd
            [pxx, psdFreq] = pwelch(lfpInState, win, noverlap, [], fs);
            pxx = log10(pxx);
            psdPow{counter}(:, isession) = pxx;
            counter = counter + 1;
        end        
    end
     
    % plot
    fh = figure;
    for istate = 1 : length(stateidx)
        subplot(1, length(stateidx), istate)
        plot(psdFreq, movmean(psdPow{istate}, 1))
        xlim([1.5 100])
        ylim([0.5 4])
        set(gca, 'xscale', 'log', 'yscale', 'linear', 'box', 'off')
        legend(sessionDate)
        legend({'K 10mg/kg', 'Saline', 'K 60mg/kg'})
        title(ss.labelNames{stateidx(istate)})
        ylabel('PSD [dB/Hz]')
        xlabel('Frequency [log(Hz)]')
    end    
    if saveFig
        figpath = fullfile(mousepath, 'graphics');
        mkdir(figpath)
        figname = fullfile(figpath, ['psd']);
        export_fig(figname, '-tif', '-transparent', '-r300')
    end   
end


% mean and median across time for all sessions (one figure)
grp = [2];          % which tetrodes to plot
unitClass = 'pyr';      % plot 'int', 'pyr', or 'all'
suFlag = 1;             % plot only su or all units
frBoundries = [0 Inf];  % include only units with fr greater than
yscale = 'log';         % log or linear

% initialize
tstamps = [];
tidx = 0;
muFr = [];
pyrMean = [];
pyrMedian = [];
intMean = [];
intMedian = [];

% gather data
for isession = 1 : nsessions
    
    % vars
    assignVars(varArray, isession)
    
    % timestamps
    tstamps = [tstamps, fr.tstamps / 60 / 60 + tidx(isession)];
    if isession == 1
        tidx(isession + 1) = fr.tstamps(end) / 60 / 60;
    else
        tidx(isession + 1) = fr.tstamps(end) / 60 / 60 + tidx(isession);
    end

    % MU
    muFr = [muFr, sr.strd];
    
    % RS
    units = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'pyr');
    pyrFr{isession} = fr.strd(units, :);
    
    % FS
    units = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'int');
    intFr{isession} = fr.strd(units, :);
end


% plot
fh = figure;
subplot(2, 1, 1)
plot(muFr')
set(gca, 'xtick', [], 'box', 'off')
ylabel('Multi-unit firing rate [Hz]')
legend
xline(tidx, '--k')
subplot(2, 1, 2)
plot(tstamps, nanmean(cell2nanmat(pyrFr, 1)))
hold on
plot(tstamps, nanmean(cell2nanmat(intFr, 1)))
xlabel('Time [h]')
ylabel('Single unit firing rate [Hz]')
set(gca, 'box', 'off')
legend('RS', 'FS')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% stacked bar plot of state duration across sessions
figFlag = 1; saveFig = true;
stateidx = [1 : 6];
sessionidx = 1 : nsessions;
[nsub] = numSubplots(length(sessionidx));
if figFlag
    clear stateDur
    for isession = sessionidx
        assignVars(varArray, isession)
        [cfg_colors, cfg_names, ~] = as_loadConfig([]);
        for istate = stateidx 
            stateDur(isession, istate) = sum(ss.epLen{istate}) /...
                length(ss.labels) * 100;
        end               
    end
    
    fh = figure;
    bh = bar(stateDur, 'stacked', 'FaceColor', 'flat');
    for isession = sessionidx
        for istate = stateidx          
            bh(istate).CData(isession, :) = cfg_colors{istate};
        end
    end
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
sessionidx = 1;
[nsub] = numSubplots(length(sessionidx));
plotStyle = 'box';      % can be 'histogram' or 'box'

if figFlag
    fh = figure;
    for isession = 1 : length(sessionidx)
        
        assignVars(varArray, sessionidx(isession))
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
            boxplot(srmat, 'PlotStyle', 'traditional', 'Whisker', 1.5);
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
        
        title(sessionDate{sessionidx(isession)})
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
sessionidx = 1;
[nsub] = numSubplots(length(sessionidx));
grp = [1 : 4];          % which tetrodes to plot
unitClass = 'int';      % plot 'int', 'pyr', or 'all'
suFlag = 1;             % plot only su or all units
frBoundries = [0 Inf];  % include only units with mean fr in these boundries

if figFlag
    fh = figure;
    clear frState
    for isession = 1 : length(sessionidx)
        assignVars(varArray, sessionidx(isession))
        units = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, unitClass);
        nunits = sum(units);
        
        for istate = stateidx
            frState{isession}(:, istate) = mean(fr.states.fr{istate}(units, :), 2, 'omitnan');
        end
        subplot(nsub(1), nsub(2), isession)
        hold on
        boxplot(frState{isession}(:, stateidx), 'PlotStyle', 'traditional', 'Whisker', 1.5);
        bh = findobj(gca, 'Tag', 'Box');
        bh = flipud(bh);
        for ibox = 1 : length(bh)
            patch(get(bh(ibox), 'XData'), get(bh(ibox), 'YData'),...
                cfg_colors{stateidx(ibox)}, 'FaceAlpha', 0.5)
        end
        xticks([])
        ylabel('SU firing rate [Hz]')
        ylim([0 15])
        title(sessionDate{sessionidx(isession)})
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


% WAKE / NREM ratio for each cell across sessions
figFlag = 1; saveFig = true;
if figFlag
    fh = figure;
    for isession = 1 : nsessions
        subplot(nsub(1), nsub(2), isession)
        

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nested function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


