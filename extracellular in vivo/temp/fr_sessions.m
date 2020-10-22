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
pcond = ["tempFlag"];
% pcond = [];
% same but imposes a negative condition
ncond = ["fepsp"];
sessionlist = 'sessionList.xlsx';       % must include extension
fs = 20000;                             % can also be loaded from datInfo

basepath = 'D:\VMs\shared\lh70\lh70_201004\lh70_201004';
% dirnames = ["lh58_200831_080808";...
%     "lh58_200901_080917";...
%     "lh58_200903_080936";...
%     "lh58_200905_080948";...
%     "lh58_200906_090914"];
% clear dirnames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get directory paths
% if exist('dirnames', 'var') && isstring(dirnames)
%     % ALT 1: user input dirnames
%     dirnames = dirnames;
% elseif ischar(sessionlist) && contains(sessionlist, 'xlsx')
%     % ALT 2: get dirnames from xlsx file
%     sessionInfo = readtable(fullfile(basepath, sessionlist));
%     icol = strcmp(sessionInfo.Properties.VariableNames, colName);
%     dirnames = string(table2cell(sessionInfo(:, icol)));
%     % check dirnames meet conditions
%     clear irow iicol
%     irow = ones(length(dirnames), 1);
%     for i = 1 : length(pcond)
%         iicol(i) = find(strcmp(sessionInfo.Properties.VariableNames, char(pcond(i))));
%         irow = irow & sessionInfo{:, iicol(i)} == 1;
%     end
%     for i = 1 : length(ncond)
%         iicol(i) = find(strcmp(sessionInfo.Properties.VariableNames, char(ncond(i))));
%         irow = irow & sessionInfo{:, iicol(i)} ~= 1;
%     end
%     dirnames = dirnames(irow);
%     dirnames(strlength(dirnames) < 1) = [];
% end
% 
% nsessions = length(dirnames);
% pathPieces = regexp(dirnames(:), '_', 'split'); % assumes filename structure: animal_date_time
% if nsessions > 1
%     sessionDate = [pathPieces{:}];
%     sessionDate = sessionDate(2 : 3 : end);
% else
%     sessionDate = pathPieces(2 : 3);
% end

% load files
if forceL || ~exist('d', 'var')
    d = cell(length(dirnames), length(vars));
    for i = 1 : nsessions
        filepath = char(fullfile(basepath, dirnames(i)));
        if ~exist(filepath, 'dir')
            warning('%s does not exist, skipping...', filepath)
            continue
        end
        cd(filepath)
        
        for ii = 1 : length(vars)
            filename = dir(['*', vars{ii}, '*']);
            if length(filename) == 1
                d{i, ii} = load(filename.name);
            else
                warning('no %s file in %s, skipping', vars{ii}, filepath)
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if forceA
    for i = 1 : nsessions
        
        % file
        filepath = char(fullfile(basepath, dirnames(i)));
        cd(filepath)
        
        % params
        session = CE_sessionTemplate(pwd, 'viaGUI', false,...
            'force', true, 'saveVar', true);
        nchans = session.extracellular.nChannels;
        fs = session.extracellular.sr;
        spkgrp = session.extracellular.spikeGroups.channels;
        
        % vars
        session = d{i, 1}.session;
        cm = d{i, 2}.cell_metrics;
        spikes = d{i, 3}.spikes;
        ss = d{i, 4}.SleepState;
        fr = d{i, 5}.fr;
        datInfo = d{i, 6}.datInfo;
%         
%         rez = runKS('basepath', filepath, 'fs', fs, 'nchans', nchans,...
%             'spkgrp', spkgrp, 'saveFinal', true, 'viaGui', false,...
%             'trange', [0 Inf], 'outFormat', 'ns');
        
%         % spikes
        fixSpkAndRes('grp', [], 'fs', fs);
        spikes = loadSpikes('session', session);
        spikes = fixCEspikes('basepath', filepath, 'saveVar', false,...
            'force', true);
        cell_metrics = ProcessCellMetrics('session', session,...
            'manualAdjustMonoSyn', false, 'summaryFigures', false,...
            'debugMode', true, 'transferFilesFromClusterpath', false,...
            'submitToDatabase', false);
%             
%         % su vs mu
        mu = [];
        spikes = cluVal('spikes', spikes, 'basepath', filepath, 'saveVar', true,...
            'saveFig', false, 'force', true, 'mu', mu, 'graphics', false,...
            'vis', 'on', 'spkgrp', spkgrp);

        % firing rate
        binsize = 60;
        winBL = [10 * 60 30 * 60];
        fr = firingRate(spikes.times, 'basepath', filepath,...
            'graphics', false, 'saveFig', false,...
            'binsize', binsize, 'saveVar', true, 'smet', 'MA',...
            'winBL', winBL);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveFig = false;

close all
grp = [1 : 4];          % which tetrodes to plot
state = [2];            % [] - all; 1 - awake; 2 - NREM
FRdata = 'strd';        % plot absolute fr or normalized
unitClass = 'int';      % plot 'int', 'pyr', or 'all'
suFlag = 1;             % plot only su or all units
minfr = 0;              % include only units with fr greater than
maxfr = 100;           % include only untis with fr lower than
Y = [0 50];              % ylim
p1 = 1;                 % firing rate vs. time, one fig per session
p2 = 0;                 % mfr across sessions, one fig

for i = 1 : nsessions
    
    session = d{i, 1}.session;
    cm = d{i, 2}.cell_metrics;
    spikes = d{i, 3}.spikes;
    ss = d{i, 4}.SleepState;
    fr = d{i, 5}.fr;
    datInfo = d{i, 6}.datInfo;
    fs = session.extracellular.sr;
   
    % cell class
    pyr = strcmp(cm.putativeCellType, 'Pyramidal Cell');
    int = strcmp(cm.putativeCellType, 'Narrow Interneuron');
    
    su = ones(1, length(spikes.ts));    % override
    if isfield(spikes, 'su') && suFlag
        su = spikes.su';
        suTxt = 'SU';
    else
        suTxt = 'MU';
    end    
    % specific grp
    grpidx = zeros(1, length(spikes.shankID));
    for ii = 1 : length(grp)
        grpidx = grpidx | spikes.shankID == grp(ii);
    end
    mfrunits = fr.mfr > minfr & fr.mfr < maxfr;
    
    if strcmp(unitClass, 'pyr')
        units = pyr & su & grpidx & mfrunits';
    elseif strcmp(unitClass, 'int')
        units = int & su & grpidx & mfrunits';
    else
        units = su & grpidx & mfrunits';     % override
    end
%     if ~any(units)
%         break
%     end
    
    switch FRdata
        case 'norm'
            data = fr.norm;
            ytxt = 'norm. MFR';
        case 'strd'
            data = fr.strd;
            ytxt = 'MFR [Hz]';
    end
    tstamps = fr.tstamps / 60;
    nsamps = cumsum(datInfo.nsamps);

    % states
    states = {ss.ints.WAKEstate, ss.ints.NREMstate, ss.ints.REMstate};
    for ii = 1 : length(states)
        tStates{ii} = InIntervals(fr.tstamps, states{ii});
        frStates{ii} = mean(data(find(units), find(tStates{ii})), 2);
    end
    % mfr in selected state across sessions (mean +- std)
    if ~isempty(state)
        mfr{i} = frStates{state};
    else
        mfr{i} = mean(data(find(units), :), 2);
    end
    totalspk(i) = length(vertcat(spikes.ts{units})) / nsamps(end) * fs;
    
    % firing rate vs. time. 1 fig per session
    if p1
        figure
        plot(tstamps, (data(units, :))')
        hold on
        medata = median((data(units, :)), 1);
        % plot(tstamps, medata, 'k', 'LineWidth', 5)
        stdshade(data(units, :), 0.3, 'k', tstamps)
        for ii = 1 : length(nsamps) - 1
            plot([nsamps(ii) nsamps(ii)] / fs / 60, ylim, '--k',...
                'LineWidth', 2)
        end
        axis tight
        fill([states{2} fliplr(states{2})]' / 60, [Y(1) Y(1) Y(2) Y(2)],...
            'b', 'FaceAlpha', 0.15,  'EdgeAlpha', 0);
        ylim(Y)
        ylabel(ytxt)
        xlabel('Time [m]')
        suptitle(dirnames{i})
        if saveFig
            figname = sprintf('%s_FRvsTime', sessionDate{i})
            figname = fullfile(basepath, figname);
            % print(fh, figname, '-dpdf', '-bestfit', '-painters');
            export_fig(figname, '-tif', '-transparent', '-r300')
        end
    end
end
clr = ['bbbkkkkkkrrr'];
if p2
    fh = figure;
    mfrmat = cell2nanmat(mfr);
    %     plot([1 : nsessions], cell2nanmat(mfr), 'ko',...
    %         'LineStyle', 'none', 'MarkerFaceColor', 'k', 'MarkerSize', 4)
    %     errorbar([1 : nsessions], nanmean(mfrmat), nanstd(mfrmat, 1),...
    %         'k', 'LineWidth', 2)
    %     hold on
    %     bh = bar([1 : nsessions], nanmean(mfrmat), 'EdgeColor', 'none');
    boxplot(mfrmat, 'PlotStyle', 'traditional');
    bh = findobj(gca, 'Tag', 'Box');
    if length(bh) == length(clr)
        clrRep = histc(clr, unique(clr));
        clear alphaIdx
        for i = 1 : length(clrRep)
            alphaIdx{i} = linspace(1, 0.3, clrRep(i));
        end
        alphaIdx = [alphaIdx{:}];
        for i = 1 : length(bh)
              patch(get(bh(i), 'XData'), get(bh(i), 'YData'),...
                clr(i), 'FaceAlpha', alphaIdx(i))
        end
    end
    ylabel(ytxt)
    xlabel('Session')
    xticks(1 : nsessions)
    xticklabels(sessionDate)
    xtickangle(45)
    box off
    title(sprintf('MFR of %s %s', unitClass, suTxt))
    if saveFig
        figname = sprintf('mfr_%s_%s', unitClass, suTxt);
        figname = fullfile(basepath, figname);
        % print(fh, figname, '-dpdf', '-bestfit', '-painters');
        export_fig(figname, '-tif', '-transparent', '-r300')
    end
end

% number of SU and MU (top - pyr; bottom - int) from selected tetrodes
p2 = 0;
grp = [1 : 8];
if p2
    clear units
    for i = 1 : nsessions
        cm = d{i, 2}.cell_metrics;
        spikes = d{i, 3}.spikes;
        pyr = strcmp(cm.putativeCellType, 'Pyramidal Cell');
        int = strcmp(cm.putativeCellType, 'Narrow Interneuron');
        grpidx = zeros(1, length(spikes.shankID));
        for ii = 1 : length(grp)
            grpidx = grpidx | spikes.shankID == grp(ii);
        end
        units(i, 1) = sum(~spikes.su' & pyr & grpidx);
        units(i, 2) = sum(spikes.su' & pyr & grpidx);
        units(i, 3) = sum(~spikes.su' & int & grpidx);
        units(i, 4) = sum(spikes.su' & int & grpidx);
        units(i, 5) = sum(grpidx) - sum(units(i, 1 : 4));
    end
    fh = figure;
    bar(units, 'stacked')
    legend({"pyr mu"; "pyr su"; "int mu"; "int su"; "others"})
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