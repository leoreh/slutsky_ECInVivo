% fr_sessions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
forceL = false;
forceA = false;

% full path and name to xls file with session metadata
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
mname = 'lh81';

% column name in xls sheet where dirnames exist
colName = 'Session';                    

% string array of variables to load
vars = ["session.mat";...
    "cell_metrics.cellinfo";...
    "spikes.cellinfo";...
    "SleepState.states";....
    "fr.mat";...
    "datInfo";...
    "AccuSleep_states"];

% column name of logical values for each session. only if true than session
% will be loaded. can be a string array and than all conditions must be
% met.
pcond = ["mancur"; "tempflag"];

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
    for i = 1 : nsessions
        
        % file
        basepath = char(fullfile(mousepath, dirnames{i}));
        basepaths{i} = fullfile(mousepath, dirnames{i});
        cd(basepath)
        [~, basename] = fileparts(basepath);
        
        % params
        session = CE_sessionTemplate(pwd, 'viaGUI', false,...
            'force', true, 'saveVar', true);
        nchans = session.extracellular.nChannels;
        fs = session.extracellular.sr;
        spkgrp = session.extracellular.spikeGroups.channels;
        
        % vars
%         session = varArray{i, 1}.session;
%         cm = varArray{i, 2}.cell_metrics;
        spikes = varArray{i, 3}.spikes;
%         if ~isempty(varArray{i, 4})
%             ss = varArray{i, 4}.SleepState;
%         end
%         fr = varArray{i, 5}.fr;
%         datInfo = varArray{i, 6}.datInfo;
%         
%         % spikes
%         fixSpkAndRes('grp', [], 'fs', fs, 'stdFactor', 0, 'nchans', nchans);
%         spikes = loadSpikes('session', session);
%         spikes = fixCEspikes('basepath', basepath, 'saveVar', false,...
%             'force', true);
%         
        spikes = cluVal('spikes', spikes, 'basepath', basepath, 'saveVar', true,...
            'saveFig', false, 'force', true, 'mu', [], 'graphics', false,...
            'vis', 'on', 'spkgrp', spkgrp);
        
cell_metrics = ProcessCellMetrics('session', session,...
    'manualAdjustMonoSyn', false, 'summaryFigures', false,...
    'debugMode', true, 'transferFilesFromClusterpath', false,...
    'submitToDatabase', false, 'getWaveformsFromDat', true);
%         cell_metrics = CellExplorer('basepath', basepath);
%         
        cc = cellclass('basepath', basepath,...
            'waves', cat(1, spikes.rawWaveform{:})', 'saveVar', true,...
            'graphics', false, 'fs', fs);
        
        % firing rate
        load([basename '.spikes.cellinfo.mat'])
        binsize = 60;
        winBL = [1 * 60 110 * 60];
        fr = firingRate(spikes.times, 'basepath', basepath,...
            'graphics', false, 'saveFig', false,...
            'binsize', binsize, 'saveVar', true, 'smet', 'MA',...
            'winBL', winBL);
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

close all
grp = [1 : 4];          % which tetrodes to plot
state = [];             % [] - all; 1 - awake; 2 - NREM
FRdata = 'strd';        % plot absolute fr or normalized
unitClass = 'pyr';      % plot 'int', 'pyr', or 'all'
suFlag = 1;             % plot only su or all units
minfr = 0;              % include only units with fr greater than
maxfr = 3000;           % include only units with fr lower than
Y = [0 10];             % ylim
yscale = 'log';         % log or linear
p1 = 1;                 % firing rate vs. time, one fig per session
p2 = 0;                 % mfr across sessions, one fig
p3 = 0;                 % firing rate vs. time, one fig for all sessions. not rubust
p4 = 0;                 % number of cells per session, one fig
plotStyle = 'box';      % for p2. can by 'bar' or 'box'
clr = ['bbbkkkkrrr'];
[nsub] = numSubplots(nsessions);

for i = 1 : nsessions
    
    % vars
    session = varArray{i, 1}.session;
    cm = varArray{i, 2}.cell_metrics;
    spikes = varArray{i, 3}.spikes;
    if ~isempty(varArray{i, 7})
        ss = varArray{i, 7}.ss;
    elseif ~isempty(varArray{i, 4})
        ss = varArray{i, 4}.SleepState;
    else
        ss = [];
    end
    fr = varArray{i, 5}.fr;
    datInfo = varArray{i, 6}.datInfo;
    fs = session.extracellular.sr;   
    
    % su vs mu
    su = ones(length(spikes.ts), 1);    % override
    if isfield(spikes, 'su') && suFlag
        su = spikes.su';
        suTxt = 'SU';
    else
        suTxt = 'MU';
    end
    
    % tetrode
    grpidx = zeros(1, length(spikes.shankID));
    for ii = 1 : length(grp)
        grpidx = grpidx | spikes.shankID == grp(ii);
    end
    mfrunits = fr.mfr > minfr & fr.mfr < maxfr;
    
    % cell class
    pyr = strcmp(cm.putativeCellType, 'Pyramidal Cell');
    wide = strcmp(cm.putativeCellType, 'Wide Interneuron');
    int = strcmp(cm.putativeCellType, 'Narrow Interneuron');
    if strcmp(unitClass, 'pyr')
        units{i} = pyr & su' & grpidx & mfrunits';
    elseif strcmp(unitClass, 'int')
        units{i} = int & su' & grpidx & mfrunits';
    else
        units{i} = su' & grpidx & mfrunits';     % override
    end
    if ~any(units{i})
        error('no units match criteria')
    end
    
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
        data = fr.states.fr{state}(units{i}, :);
        tstamps = fr.states.tstamps{state};
        txt = sprintf('mFR of %s %s in %s', unitClass, suTxt, fr.states.statenames{state});
    else
        data = fr.(FRdata)(units{i}, :);
        tstamps = fr.tstamps;
        txt = sprintf('mFR of %s %s', unitClass, suTxt);
    end
    mfr{i} = mean(data, 2);
    
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
        suptitle(dirnames{i})
        if saveFig
            if nsessions == 1
                clear sessionDate
                sessionDate{1} = pathPieces{2};
            end
            figname = sprintf('%s_%s_FRvsTime', sessionDate{i}, unitClass)
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
    
    % fr vs. time; one figure for all sessoins
    k = [1, 2, 3, 4]';
    %     titletxt = ["lh69_200907_morning"; "lh69_200908_morning";...
    %          "lh69_200908_evening"];
    if p3
        subplot(2, 2, k(i))
        plot(tstamps / 60, mean(data(units{i}, :)), 'k', 'LineWidth', 2)
        hold on
        if i == 3
            for ii = 2
                plot([nsamps(ii) nsamps(ii)] / fs / 60 / 60, Y, '--k',...
                    'LineWidth', 2)
            end
        end
        axis tight
        fill([states{2} fliplr(states{2})]' / 60 / 60, [Y(1) Y(1) Y(2) Y(2)],...
            'b', 'FaceAlpha', 0.15,  'EdgeAlpha', 0);
        ylim(Y)
        ylabel(ytxt)
        xlabel('Time [h]')
        title(dirnames{i})
        if saveFig
            figname = sprintf('LTP')
            figname = fullfile(mousepath, figname);
            % print(fh, figname, '-dpdf', '-bestfit', '-painters');
            export_fig(figname, '-tif', '-transparent', '-r300')
        end
    end
end

if p2
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
                for i = 1 : length(clrRep)
                    alphaIdx{i} = linspace(0.8, 0.3, clrRep(i));
                end
                alphaIdx = [alphaIdx{:}];
                for i = 1 : length(bh)
                    patch(get(bh(i), 'XData'), get(bh(i), 'YData'),...
                        clr(i), 'FaceAlpha', alphaIdx(i))
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
if p4
    clear units
    for i = 1 : nsessions
        cm = varArray{i, 2}.cell_metrics;
        spikes = varArray{i, 3}.spikes;
        pyr = strcmp(cm.putativeCellType, 'Pyramidal Cell');
        wide = strcmp(cm.putativeCellType, 'Wide Interneuron');
        pyr = pyr | wide;
        int = strcmp(cm.putativeCellType, 'Narrow Interneuron');
        grpidx = zeros(1, length(spikes.shankID));
        for ii = 1 : length(grp)
            grpidx = grpidx | spikes.shankID == grp(ii);
        end
        units(i, 1) = sum(~spikes.su & pyr & grpidx);
        units(i, 2) = sum(spikes.su & pyr & grpidx);
        units(i, 3) = sum(~spikes.su & int & grpidx);
        units(i, 4) = sum(spikes.su & int & grpidx);
        %         units(i, 5) = sum(grpidx) - sum(units(i, 1 : 4));
    end
    fh = figure;
    bar(units, 'stacked')
    hold on    
    if nsessions == length(clr)
        idx = [find(clr == 'k', 1) find(clr == 'r', 1)] - 0.5;
        yLimit = ylim;
        patch([idx idx(2) idx(1)],...
            [yLimit(1) yLimit(1) yLimit(2) yLimit(2)],...
            'b', 'FaceAlpha', 0.1)
    end
    
    legend({"pyr mu"; "pyr su"; "int mu"; "int su"; "CNO"})
    xticks(1 : nsessions)
    xticklabels(sessionDate)
    title('Number of Units')
    xlabel('Session')
    ylabel('No. units')
    box off
    
    if saveFig
        figname = fullfile(mousepath, 'UnitsDetected');
        % print(fh, figname, '-dpdf', '-bestfit', '-painters');
        export_fig(figname, '-tif', '-transparent', '-r300')
    end
end

% mean firing rate in states per unit
p3 = 0;
k = 1;
% data = [];
stateidx = 1 : 4;
if p3
    fh = figure;
    for i = 1 : nsessions
           % vars
           session = varArray{i, 1}.session;
           cm = varArray{i, 2}.cell_metrics;
           spikes = varArray{i, 3}.spikes;
           fr = varArray{i, 5}.fr;
           datInfo = varArray{i, 6}.datInfo;
           fs = session.extracellular.sr;
           ss = varArray{i, 7}.ss;
           nunits = length(spikes.ts);
           
           clear data
           for ii = 1 : length(ss.labelNames)
               data(:, ii) = mean(fr.states.fr{ii}(units{i}, :), 2);
%                if ~isempty(data)
%                    h = histogram(data, round(nunits / 4), 'EdgeAlpha', 0,...
%                        'FaceAlpha', 0.3, 'Normalization', 'count');
%                end
           end
           
           subplot(nsub(1), nsub(2), k)
           plot(data(:, stateidx)')
           hold on
           plot(mean(data(:, stateidx)), 'k', 'LineWidth', 2)
           xticks([1 : length(stateidx)])
           xticklabels(fr.states.stateNames(stateidx))
        ylabel('Mean Firing Rate [Hz]')
        xlabel('State')
        if k == 1
            legend(fr.states.stateNames(stateidx))
        end
        
%         xlim([0 150])
        if i == 1
            legend(ss.labelNames)
        end
        title(dirnames{i})
        k = k + 1;
    end
end