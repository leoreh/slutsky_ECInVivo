% sr_sessions
% spike rate per tetrode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
forceL = false;
forceA = false;
basepath = 'D:\VMs\shared\lh70';

% should allow user to input varName or columnn index
colName = 'Session';                    % column name in xls sheet where dirnames exist
% string array of variables to load
vars = ["session.mat";...
    "SleepState.states";....
    "datInfo";
    ".sr.mat"];
% column name of logical values for each session. only if true than session
% will be loaded. can be a string array and than all conditions must be
% met.
pcond = ["tempflag"];
% pcond = [];
% same but imposes a negative condition
ncond = ["fepsp"];
sessionlist = 'sessionList.xlsx';       % must include extension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('varArray', 'var') && ~forceL
    [varArray, dirnames] = getSessionVars('basepath', basepath, 'vars', vars,...
        'pcond', pcond, 'ncond', ncond, 'sortDir', false, 'dirnames', []);
end
nsessions = length(dirnames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if forceA
    for i = 1 : nsessions
        
        % file
        filepath = char(fullfile(basepath, dirnames{i}));
        [~, basename] = fileparts(filepath);
        cd(filepath)
        
        % params
        session = CE_sessionTemplate(pwd, 'viaGUI', false,...
            'force', true, 'saveVar', true);
        nchans = session.extracellular.nChannels;
        fs = session.extracellular.sr;
        spkgrp = session.extracellular.spikeGroups.channels;
                
        % detect spikes
        spktimesWh('basepath', filepath, 'fs', fs, 'nchans', nchans,...
            'spkgrp', spkgrp, 'saveVar', true, 'saveWh', false,...
            'graphics', false)
        
        % create ns files for sorting
        spktimes2ks('basepath', filepath, 'fs', fs,...
            'nchans', nchans, 'spkgrp', spkgrp, 'mkClu', true,...
            'dur', 240, 't', '0200', 'psamp', [], 'grps', [1 : length(spkgrp)]);
        
        % firing rate per tetrode. note that using times2rate requires special care
        % becasue spktimes is given in samples and not seconds
        binsize = 60 * fs;
        winCalc = [0 Inf];
        load([basename '.spktimes.mat'])
        [sr.strd, sr.edges, sr.tstamps] = times2rate(spktimes, 'binsize', binsize,...
            'winCalc', winCalc, 'c2r', false);
        % convert counts to rate
        sr.strd = sr.strd ./ (diff(sr.edges) / fs);
        % fix tstamps
        sr.tstamps = sr.tstamps / binsize;        
        save(fullfile(filepath, [basename '.sr.mat']), 'sr')
    end
end

% params
session = varArray{1, 1}.session;
% session = CE_sessionTemplate(pwd, 'viaGUI', false,...
%     'force', true, 'saveVar', true);
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
 
% close all
grp = [1 : 4];          % which tetrodes to plot
state = [2];            % [] - all; 1 - awake; 2 - NREM
Y = [0 200];             % ylim
p1 = 1;                 % firing rate vs. time, one fig per session
p2 = 0;                 % mfr across sessions, one fig
p3 = 0;                 % firing rate vs. time, one fig for all sessions. not rubust 

for i = 1 : nsessions
    
    session = varArray{i, 1}.session;
    ss = varArray{i, 2}.SleepState;
    datInfo = varArray{i, 3}.datInfo;
    sr = varArray{i, 4}.sr;
   
    data = sr.strd;
    tstamps = sr.tstamps;
    if isfield(datInfo, 'nsamps')
        nsamps = cumsum(datInfo.nsamps);
    else
        nsamps = 120 * 60 * fs;
    end
    
    % states
    states = {ss.ints.WAKEstate, ss.ints.NREMstate, ss.ints.REMstate};
    for ii = 1 : length(states)
        tStates{ii} = InIntervals(sr.tstamps * 60, states{ii});
        frStates{ii} = mean(data(grp, find(tStates{ii})), 2);
    end
    % mfr in selected state across sessions (mean +- std)
    if ~isempty(state)
        mfr{i} = frStates{state};
    else
        mfr{i} = mean(data(grp, :), 2);
    end
    
    % firing rate vs. time. 1 fig per session
    if p1
        figure
        plot(tstamps, (data(grp, :))')
        hold on
        % medata = median((data(grp, :)), 1);
        % plot(tstamps, medata, 'k', 'LineWidth', 5)
        stdshade(data(grp, :), 0.3, 'k', tstamps)
        for ii = 1 : length(nsamps) - 1
            plot([nsamps(ii) nsamps(ii)] / fs / 60, ylim, '--k',...
                'LineWidth', 2)
        end
        axis tight
        fill([states{2} fliplr(states{2})]' / 60, [Y(1) Y(1) Y(2) Y(2)],...
            'b', 'FaceAlpha', 0.15,  'EdgeAlpha', 0);
        ylim(Y)
        ylabel('Spike Rate [Hz]')
        xlabel('Time [m]')
        suptitle(dirnames{i})
        if saveFig
            figname = sprintf('%s_FRvsTime', sessionDate{i})
            figname = fullfile(basepath, figname);
            % print(fh, figname, '-dpdf', '-bestfit', '-painters');
            export_fig(figname, '-tif', '-transparent', '-r300')
        end
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
%     boxplot(mfrmat, 'PlotStyle', 'traditional');
    bar(nanmean(mfrmat))
    hold on
    errorbar(1 : nsessions, nanmean(mfrmat), nanstd(mfrmat))
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
    ylabel('Spike Rate [Hz]')
    xlabel('Session')
    xticks(1 : nsessions)
    xticklabels(sessionDate)
    xtickangle(45)
    box off
    title(sprintf('MSR of Tetrode #%s', num2str(grp)))
    if saveFig
        figname = sprintf('MSR of Tetrode #%s', num2str(grp));
        figname = fullfile(basepath, figname);
        % print(fh, figname, '-dpdf', '-bestfit', '-painters');
        export_fig(figname, '-tif', '-transparent', '-r300')
    end
end

% number of SU and MU (top - pyr; bottom - int) from selected tetrodes
p2 = 0;
% grp = [1 : 8];
if p2
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
        units(i, 1) = sum(~spikes.su' & pyr & grpidx);
        units(i, 2) = sum(spikes.su' & pyr & grpidx);
        units(i, 3) = sum(~spikes.su' & int & grpidx);
        units(i, 4) = sum(spikes.su' & int & grpidx);
        %         units(i, 5) = sum(grpidx) - sum(units(i, 1 : 4));
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