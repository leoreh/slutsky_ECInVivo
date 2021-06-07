% fr_sessions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
forceL = false;
forceA = false;

% full path and name to xls file with session metadata
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
mname = 'lh91';

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
pcond = ["tempflag"];

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
        %
        %         cc = cellclass('basepath', basepath,...
        %             'waves', cat(1, spikes.rawWaveform{:})', 'saveVar', true,...
        %             'graphics', false, 'fs', fs);
        
        %         spikes = cluVal('spikes', spikes, 'basepath', basepath, 'saveVar', true,...
        %             'saveFig', false, 'force', true, 'mu', [], 'graphics', false,...
        %             'vis', 'on', 'spkgrp', spkgrp);
        
        % firing rate
        load([basename '.spikes.cellinfo.mat'])
        binsize = 60;
        winBL = [1 * 60 110 * 60];
        fr = firingRate(spikes.times, 'basepath', basepath,...
            'graphics', false, 'saveFig', false,...
            'binsize', binsize, 'saveVar', true, 'smet', 'MA',...
            'winBL', winBL);
        
        % spike rate per tetrode. note that using firingRate requires
        % special care becasue spktimes is given in samples and not seconds
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
saveFig = false;
pathPieces = regexp(dirnames(:), '_', 'split'); % assumes filename structure: animal_date_time
sessionDate = [pathPieces{:}];
sessionDate = sessionDate(2 : 3 : end);
setMatlabGraphics(false)
[nsub] = numSubplots(nsessions);

% close all
grp = [1 : 4];          % which tetrodes to plot
state = [];             % [] - all; 1 - awake; 2 - NREM
FRdata = 'strd';        % plot absolute fr or normalized
unitClass = 'int';      % plot 'int', 'pyr', or 'all'
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

% mean firing rate in states per unit
figFlag = 0;
stateidx = [1, 4];
unitClass = 'pyr';
clear stateFrUnits
if figFlag
    fh = figure;
    for isession = 1 : nsessions
        
        assignVars(varArray, isession)
        units = selectUnits(spikes, cm, fr, 1, grp, frBoundries, unitClass);
        nunits = length(units);
        
        subplot(2, nsessions, isession)
        hold on
        k = 1;
        for istate = stateidx
            % select specific times
            idx = fr.states.tstamps{istate} > 2 * 60 * 60 &...
                fr.states.tstamps{istate} < 6 * 60 * 60;
            frmat = fr.states.fr{istate}(units, idx);
            stateFrUnits{isession}(:, k) = mean(frmat, 2);
            k = k + 1;
        end
        
        %         for istate = stateidx
        %             if ~isempty(stateFrUnits{isession})
        %             h = histogram(stateFrUnits{isession}, round(nunits / 4), 'EdgeAlpha', 0,...
        %                 'FaceAlpha', 0.3, 'Normalization', 'countdensity');
        %             end
        %         end
        %         set(gca, 'xscale', 'log', 'box', 'off')
        %         xlabel('Firing Rate [Hz]')
        %         ylabel('No. Units')
        %         xlim([0 100])
        %         if isession == 1
        %             legend(ss.labelNames(stateidx))
        %         end
        
        plot(stateFrUnits{isession}')
        hold on
        plot(median(stateFrUnits{isession}), 'k', 'LineWidth', 2)
        xLimit = xlim;
        xLimit = [xLimit(1) - 0.5, xLimit(2) + 0.5];
        xlim(xLimit)
        xticks([1 : length(stateidx)])
        xticklabels(fr.states.stateNames(stateidx))
        ylabel('Mean Firing Rate [Hz]')
        xlabel('State')
        set(gca, 'yscale', 'log', 'box', 'off')
        ylim([0.1 10])      
        title(sessionDate{isession})
        k = k + 1;
    end
    
    % percent change across sessions
    subplot(2, nsessions, nsessions + 1 : nsessions * 2)
    clear percentChange
    for isession = 1 : nsessions
        percentChange{isession} = stateFrUnits{isession}(:, :) ./ stateFrUnits{isession}(:, 1);
        percentChange{isession} = percentChange{isession}(:, 2 : end);
    end
    data = cell2nanmat(percentChange);
    boxplot(data, 'PlotStyle', 'traditional');
    ylabel('FR ratio NREM / WAKE')
    
    if saveFig
        figpath = fullfile(mousepath, 'graphics');
        mkdir(figpath)
        figname = fullfile(figpath, ['states_', unitClass]);
        export_fig(figname, '-tif', '-transparent', '-r300')
    end
end

% histogram of bins for each state for multiunit activity 
figFlag = 0;
grp = 1 : 4;
stateidx = [1, 4];
if figFlag
    fh = figure;
    for isession = 1 : nsessions
        
        assignVars(varArray, isession)
        units = selectUnits(spikes, cm, fr, 1, grp, frBoundries, unitClass);
        nunits = length(units);
        
        subplot(nsub(1), nsub(2), isession)
        hold on
        for istate = stateidx
            tidx = sr.states.tstamps{istate} > 2 * 60 * 60;
            srmat = sr.states.fr{istate}(:, tidx);
            nbins = min([length(sr.states.fr{istate}), 250]);
            binidx = randperm(length(sr.states.fr{istate}), nbins);
            binstates{isession, istate} = mean(sr.states.fr{istate}(grp, binidx), 1);
            if ~isempty(binstates{isession, istate}) && nbins > 10
                h = histogram(binstates{isession, istate}, round(nbins / 10), 'EdgeAlpha', 0,...
                    'FaceAlpha', 0.3, 'Normalization', 'count');
            end
        end
        
        set(gca, 'xscale', 'log')
        xlim([10 250])
        ylim([0 60])
        ylabel('Counts')
        xlabel('Spike Rate [Hz]')
        title(sessionDate{isession})   
        if isession == 1
            legend(ss.labelNames(stateidx))
        end
    end
    
    if saveFig
        figpath = fullfile(mousepath, 'graphics');
        mkdir(figpath)
        figname = fullfile(figpath, ['states_mu']);
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
grp = [1 : 4];          % which tetrodes to plot
unitClass = 'int';      % plot 'int', 'pyr', or 'all'
suFlag = 1;             % plot only su or all units
frBoundries = [0 Inf];  % include only units with fr greater than
yscale = 'log';         % log or linear

% initialize
tstamps = [];
tidx = [];
muFr = [];
pyrMean = [];
pyrMedian = [];
intMean = [];
intMedian = [];

% gather data
% for i = 1 : nsessions
%     
%     % vars
%     assignVars(varArray, isession)
%     
%     % timestamps
%     tstamps = [tstamps, fr.tstamps / 60 / 60 + tidx(isession)];
%     if isession == 1
%         tidx(isession) = fr.tstamps(end) / 60 / 60;
%     else
%         tidx(isession) = fr.tstamps(end) / 60 / 60 + tidx(isession - 1);
%     end
%     
%     % MU
%     muFr = [muFr, sr.strd];
%     
%     % RS
%     units = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'pyr');
%     pyrMean = [pyrMean, mean(fr.strd(units, :), 1)];
%     pyrMedian = [pyrMedian, median(fr.strd(units, :), 1)];
%     
%     % FS
%     units = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'int');
%     intMean = [intMean, mean(fr.strd(units, :), 1)];
%     intMedian = [intMedian, median(fr.strd(units, :), 1)];
%     
% end

% plot
% 
% fh = figure;
% subplot(2, 1, 1)
% plot(muFr')
% set(gca, 'xtick', [], 'box', 'off')
% ylabel('Multi-unit firing rate [Hz]')
% subplot(2, 1, 2)
% plot(tstamps, pyrMean)
% hold on
% plot(tstamps, intMean)
% xlabel('Time [h]')
% ylabel('Single unit firing rate [Hz]')
% set(gca, 'box', 'off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nested function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


