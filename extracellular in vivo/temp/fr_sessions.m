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
ncond = ["fix"];
ncond = ["fepsp"];
sessionlist = 'sessionList.xlsx';       % must include extension
fs = 20000;                             % can also be loaded from datInfo

basepath = 'G:\Data\Processed\lh58';
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
if exist('dirnames', 'var') && isstring(dirnames)
    % ALT 1: user input dirnames
    dirnames = dirnames;
elseif ischar(sessionlist) && contains(sessionlist, 'xlsx')
    % ALT 2: get dirnames from xlsx file
    sessionInfo = readtable(fullfile(basepath, sessionlist));
    icol = strcmp(sessionInfo.Properties.VariableNames, colName);
    dirnames = string(table2cell(sessionInfo(:, icol)));
    % check dirnames meet conditions
    clear irow iicol
    irow = ones(length(dirnames), 1);
    for i = 1 : length(pcond)
        iicol(i) = find(strcmp(sessionInfo.Properties.VariableNames, char(pcond(i))));
        irow = irow & sessionInfo{:, iicol(i)} == 1;
    end
    for i = 1 : length(ncond)
        iicol(i) = find(strcmp(sessionInfo.Properties.VariableNames, char(ncond(i))));
        irow = irow & sessionInfo{:, iicol(i)} ~= 1;
    end
    dirnames = dirnames(irow);
    dirnames(strlength(dirnames) < 1) = [];
end

nsessions = length(dirnames);
pathPieces = regexp(dirnames(:), '_', 'split'); % assumes filename structure: animal_date_time
sessionDate = [pathPieces{:}];
sessionDate = sessionDate(2 : 3 : end);

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
        
        % vars
        session = d{i, 1}.session;
        cm = d{i, 2}.cell_metrics;
        spikes = d{i, 3}.spikes;
        ss = d{i, 4}.SleepState;
        fr = d{i, 5}.fr;
        datInfo = d{i, 6}.datInfo;
        
        % params
        nchans = session.extracellular.nChannels;
        fs = session.extracellular.sr;
        spkgrp = session.extracellular.spikeGroups.channels;
        
        % spikes
%         fixSpkAndRes('grp', [], 'fs', fs);
%         spikes = loadSpikes('session', session);
%         spikes = fixCEspikes('basepath', filepath, 'saveVar', false,...
%             'force', true);
%         cell_metrics = ProcessCellMetrics('session', session,...
%             'manualAdjustMonoSyn', false, 'summaryFigures', false,...
%             'debugMode', true, 'transferFilesFromClusterpath', false,...
%             'submitToDatabase', false);
            
        % su vs mu
        mu = [];
        spikes = cluVal('spikes', spikes, 'basepath', filepath, 'saveVar', true,...
            'saveFig', false, 'force', true, 'mu', mu, 'graphics', false,...
            'vis', 'on', 'spkgrp', spkgrp);

        % firing rate
%         binsize = 300;
%         winBL = [10 * 60 30 * 60];
%         fr = firingRate(spikes.times, 'basepath', filepath,...
%             'graphics', false, 'saveFig', false,...
%             'binsize', binsize, 'saveVar', true, 'smet', 'MA',...
%             'winBL', winBL);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
grp = [1 : 8];          % which tetrodes to plot
state = [];              % [] - all; 1 - awake; 2 - NREM
FRdata = 'strd';        % plot absolute fr or normalized
unitClass = 'pyr';      % plot 'int', 'pyr', or 'all'
suFlag = true;         % plot only su or all units
minfr = 0;              % include only units with fr greater than
p1 = false;             % firing rate vs. time, one fig per session
p2 = true;              % mfr across sessions, one fig

for i = 1 : nsessions
    
    session = d{i, 1}.session;
    cm = d{i, 2}.cell_metrics;
    spikes = d{i, 3}.spikes;
    ss = d{i, 4}.SleepState;
    fr = d{i, 5}.fr;
    datInfo = d{i, 6}.datInfo;
       
    % cell class
    pyr = strcmp(cm.putativeCellType, 'Pyramidal Cell');
    int = strcmp(cm.putativeCellType, 'Narrow Interneuron');
    
    su = ones(1, length(spikes.ts));    % override
    if isfield(spikes, 'su') && suFlag
        su = spikes.su';
    end    
    % specific grp
    grpidx = zeros(1, length(spikes.shankID));
    for ii = 1 : length(grp)
        grpidx = grpidx | spikes.shankID == grp(ii);
    end
    mfrunits = fr.mfr > minfr;
    
    if strcmp(unitClass, 'pyr')
        units = pyr & su & grpidx & mfrunits';
    elseif strcmp(unitClass, 'int')
        units = int & su & grpidx & mfrunits';
    else
        units = su & grpidx & mfrunits';     % override
    end
    
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
    
    % firing rate vs. time. 1 fig per session
    if p1
        figure
        plot(tstamps, (data(units, :))')
        hold on
        medata = median((data(units, :)), 1);
        % plot(tstamps, medata, 'k', 'LineWidth', 5)
        stdshade(data(units, :), 0.3, 'k', tstamps)
        for ii = 1 : length(nsamps) - 1
            plot([nsamps(ii) nsamps(ii)] / fs / 60, ylim, '--k')
        end
        axis tight
        Y = ylim;
        fill([states{2} fliplr(states{2})]' / 60, [Y(1) Y(1) Y(2) Y(2)],...
            'b', 'FaceAlpha', 0.15,  'EdgeAlpha', 0);
        ylim([0 3])
        ylabel(ytxt)
        suptitle(dirnames{i})
    end
    
    % mfr in selected state across sessions (mean +- std)
    if ~isempty(state)
        mfr{i} = mean(frStates{state});
        sfr(i) = std(frStates{state});
    else
        mfr{i} = mean(data(find(units), :), 2);
        sfr{i} = std(data(find(units), :), [], 2);
    end
end
if p2
    if i == 1
        fh = figure;
    end
    plot([1 : nsessions], cell2nanmat(mfr), 'ko', 'LineStyle', 'none',)
    %     bar([1 : nsessions], mfr)
%     hold on
%     errorbar([1 : nsessions], mfr, sfr, 'vertical')
%     xticks(1 : nsessions)
%     xticklabels(sessionDate)
%     box off
%     title(sprintf('MFR of %s units', unitClass))
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
    figure
    bar(units, 'stacked')
    legend({"pyr mu"; "pyr su"; "int mu"; "int su"; "others"})
    xticks(1 : nsessions)
    xticklabels(sessionDate)
    title('Number of Units')
    xlabel('Session')
    ylabel('No. units')
    box off
end