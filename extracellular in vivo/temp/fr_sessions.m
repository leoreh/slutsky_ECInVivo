% fr_sessions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
forceL = false;
forceA = false;

% should allow user to input varName or columnn index
colName = 'Session';                    % column name in xls sheet where dirnames exist
% string array of variables to load
vars = ["session.mat";...
    "cell_metrics.cellinfo";...
    "spikes.cellinfo";...
    "SleepState.states";....
    "fr.mat";...
    "datInfo"];
% column name of logical values for each session. only if true than session
% will be loaded. can be a string array and than all conditions must be
% met.
pcond = ["tempflag"; "mancur"];
pcond = "";
% same but imposes a negative condition
ncond = ["fepsp"];
sessionlist = 'sessionList.xlsx';       % must include extension

basepath = 'D:\VMs\shared\lh70';

dirnames = [];
dirnames = ["lh70_201008_0825"];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('varArray', 'var') && ~forceL
    [varArray, dirnames] = getSessionVars('basepath', basepath, 'vars', vars,...
        'pcond', pcond, 'ncond', ncond, 'sortDir', false, 'dirnames', dirnames);
end
nsessions = length(dirnames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if forceA
    for i = 1 : nsessions
        
        % file
        filepath = char(fullfile(basepath, dirnames{i}));
        basepaths{i} = fullfile(basepath, dirnames{i});
        cd(filepath)
        
        % params
        session = CE_sessionTemplate(pwd, 'viaGUI', false,...
            'force', true, 'saveVar', true);
        nchans = session.extracellular.nChannels;
        fs = session.extracellular.sr;
        spkgrp = session.extracellular.spikeGroups.channels;
        
        % vars
        %         session = varArray{i, 1}.session;
        %         cm = varArray{i, 2}.cell_metrics;
        %         spikes = varArray{i, 3}.spikes;
        %         ss = varArray{i, 4}.SleepState;
        %         fr = varArray{i, 5}.fr;
        %         datInfo = varArray{i, 6}.datInfo;
        
        % spikes
        fixSpkAndRes('grp', [1], 'fs', fs);
        spikes = loadSpikes('session', session);
        spikes = fixCEspikes('basepath', filepath, 'saveVar', false,...
            'force', true);
        
        spikes = cluVal('spikes', spikes, 'basepath', filepath, 'saveVar', true,...
            'saveFig', false, 'force', true, 'mu', [], 'graphics', false,...
            'vis', 'on', 'spkgrp', spkgrp);
        
        cell_metrics = ProcessCellMetrics('session', session,...
            'manualAdjustMonoSyn', false, 'summaryFigures', false,...
            'debugMode', true, 'transferFilesFromClusterpath', false,...
            'submitToDatabase', false, 'spikes', spikes);
        % cell_metrics = CellExplorer('basepath', filepath);
        
        cc = cellclass('basepath', filepath,...
            'waves', cat(1, spikes.rawWaveform{:})', 'saveVar', true,...
            'graphics', false);
        
        % firing rate
        binsize = 60;
        winBL = [10 * 60 30 * 60];
        fr = firingRate(spikes.times, 'basepath', filepath,...
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
state = [2];             % [] - all; 1 - awake; 2 - NREM
FRdata = 'strd';        % plot absolute fr or normalized
unitClass = 'pyr';      % plot 'int', 'pyr', or 'all'
suFlag = 1;             % plot only su or all units
minfr = 0;            % include only units with fr greater than
maxfr = 30;              % include only untis with fr lower than
Y = [0 25];             % ylim
p1 = 1;                 % firing rate vs. time, one fig per session
p2 = 0;                 % mfr across sessions, one fig
p3 = 0;                 % firing rate vs. time, one fig for all sessions. not rubust
p4 = 0;                 % number of cells per session, one fig
plotStyle = 'box';      % for p2. can by 'bar' or 'box'
clr = ['bbbkkkkrrr'];

for i = 1 : nsessions
    
    % vars
    session = varArray{i, 1}.session;
    cm = varArray{i, 2}.cell_metrics;
    spikes = varArray{i, 3}.spikes;
    ss = varArray{i, 4}.SleepState;
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
        units = pyr & su' & grpidx & mfrunits';
    elseif strcmp(unitClass, 'int')
        units = int & su' & grpidx & mfrunits';
    else
        units = su' & grpidx & mfrunits';     % override
    end
    if ~any(units)
        error('no units match criteria')
    end
    
    % states
    states = {ss.ints.WAKEstate, ss.ints.NREMstate, ss.ints.REMstate};
    if ~isempty(state) && state > 0
        data = fr.states.fr{state}(units, :);
        tstamps = fr.states.tstamps{state};
        txt = sprintf('mFR of %s %s in %s', unitClass, suTxt, fr.states.statenames{state});
    else
        data = fr.strd(units, :);
        tstamps = fr.tstamps;
        txt = sprintf('mFR of %s %s', unitClass, suTxt);
    end
    mfr{i} = mean(data, 2);
    
    %%%
    mfr{i} = cm.burstIndex_NREMstate(units);
    %%%
    
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
        medata = mean((data), 1);
        plot(tstamps / 60, medata, 'k', 'LineWidth', 4)
%         stdshade(data, 0.3, 'k', tstamps / 60)
        for ii = 1 : length(nsamps) - 1
            plot([nsamps(ii) nsamps(ii)] / fs / 60, ylim, '--k',...
                'LineWidth', 2)
        end
        axis tight
        fill([states{2} fliplr(states{2})]' / 60, [Y(1) Y(1) Y(2) Y(2)],...
            'b', 'FaceAlpha', 0.15,  'EdgeAlpha', 0);
        ylim(Y)
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
            figname = fullfile(basepath, figname);
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
        plot(tstamps / 60, mean(data(units, :)), 'k', 'LineWidth', 2)
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
            figname = fullfile(basepath, figname);
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
            boxplot(mfrmat, 'PlotStyle', 'traditional',...
                'DataLim', [minfr maxfr]);
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
        figname = fullfile(basepath, figname);
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
        figname = fullfile(basepath, 'UnitsDetected');
        % print(fh, figname, '-dpdf', '-bestfit', '-painters');
        export_fig(figname, '-tif', '-transparent', '-r300')
    end
end