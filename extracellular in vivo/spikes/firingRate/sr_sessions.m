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
    ".sr.mat";
    "spktimes.mat"];
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
        
        % sleep states
%         badch = setdiff([session.extracellular.electrodeGroups.channels{:}],...
%             [session.extracellular.spikeGroups.channels{:}]);
%         SleepScoreMaster(filepath, 'noPrompts', true,...
%             'rejectChannels', badch, 'overwrite', true,...
%             'NotchTheta', true, 'ignoretime', badstamps)
%         TheStateEditor(fullfile(filepath, basename))
%         % find bad times from sw spect    
%         badtimes = SleepState.detectorinfo.StatePlotMaterials.swFFTspec(1, :) > 0.6e4;
%         figure, plot(SleepState.detectorinfo.StatePlotMaterials.swFFTspec(1, :))
%         badstamps = binary2epochs('vec', badtimes', 'minDur', 2, 'maxDur', 200000000,...
%             'interDur', 10, 'exclude', false);
        
        % detect spikes
%         [spktimes, ~] = spktimesWh('basepath', filepath, 'fs', fs, 'nchans', nchans,...
%             'spkgrp', spkgrp, 'saveVar', true, 'saveWh', false,...
%             'graphics', false);
        
        % create ns files for sorting
        spktimes2ks('basepath', filepath, 'fs', fs,...
            'nchans', nchans, 'spkgrp', spkgrp, 'mkClu', true,...
            'dur', 240, 't', '020000', 'psamp', [], 'grps', [1 : length(spkgrp)]);
        
        % spike rate per tetrode. note that using firingRate requires
        % special care becasue spktimes is given in samples and not seconds
        load(fullfile(filepath, [basename '.spktimes.mat']))
        for ii = 1 : length(spkgrp)
            spktimes{ii} = spktimes{ii} / fs;
        end
        binsize = 60;
        sr = firingRate(spktimes, 'basepath', filepath,...
            'graphics', false, 'saveFig', false,...
            'binsize', binsize, 'saveVar', 'sr', 'smet', 'none',...
            'winBL', [0 Inf]);     
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
grp = [1 : 4];     % which tetrodes to plot
state = [2];            % [] - all; 1 - awake; 2 - NREM; 3 - REM
Y = [0 300];            % ylim
p1 = 0;                 % firing rate vs. time, one fig per session
p2 = 1;                 % mfr across sessions, one fig
p3 = 0;                 % spike rate vs. time, one fig for all sessions
plotStyle = 'bar';      % for p2. can by 'bar' or 'box'
clr = ['bbbbbkkkkkkkkrrrrrrrr'];        % color sessions

for i = 1 : nsessions
    
    session = varArray{i, 1}.session;
    ss = varArray{i, 2}.SleepState;
    datInfo = varArray{i, 3}.datInfo;
    sr = varArray{i, 4}.fr;
   
    tstamps = sr.tstamps;
    if isfield(datInfo, 'nsamps')
        nsamps = cumsum(datInfo.nsamps);
    else
        nsamps = 120 * 60 * fs;
    end
      
    % mfr in selected state across sessions (mean +- std)
    if ~isempty(state)
        mfr{i} = mean(sr.states.fr{state}(grp, :), 2);
                
%         mfr{i} = std(sr.states.fr{state}(grp, :), [], 2);

        txt = sprintf('MSR of T#%s in %s', num2str(grp), sr.states.statenames{state});
    else
        mfr{i} = mean(sr.strd(grp, :), 2);
        
%         mfr{i} = std(sr.strd(grp, :), [], 2);

        txt = sprintf('MSR of T#%s', num2str(grp));
    end
    
    % firing rate vs. time. 1 fig per session
    if p1
        figure
        plot(sr.tstamps, (sr.strd(grp, :))')
        hold on
        % medata = median((data(grp, :)), 1);
        % plot(tstamps, medata, 'k', 'LineWidth', 5)
        stdshade(sr.strd(grp, :), 0.3, 'k', tstamps)
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
end

if p2
    fh = figure;
    mfrmat = cell2nanmat(mfr);
    switch plotStyle
        case 'bar'
            plot([1 : nsessions], mfrmat)
            hold on
            bh = bar(nanmean(mfrmat), 'FaceAlpha', 0.1, 'FaceColor', 'k');
            hold on
            % errorbar(1 : nsessions, nanmean(mfrmat), nanstd(mfrmat), 'k')
            idx = [find(clr == 'k', 1) find(clr == 'r', 1)] - 0.5;
            yLimit = ylim;
            patch([idx idx(2) idx(1)],...
                [yLimit(1) yLimit(1) yLimit(2) yLimit(2)],...
                'b', 'FaceAlpha', 0.1)
            legend(split(num2str(grp)))
        case 'box'
            boxplot(mfrmat, 'PlotStyle', 'traditional');
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

% spike rate vs. time across all sessions (1 fig)
% state = [];            
if p3  
    msr = [];
    tsr = [];
    NREMstates = [];
    prev_t = 0;
    for i = 1 : nsessions
        session = varArray{i, 1}.session;
        ss = varArray{i, 2}.SleepState;
        datInfo = varArray{i, 3}.datInfo;
        sr = varArray{i, 4}.fr;
        st = varArray{i, 5}.spktimes;
        filepath = char(fullfile(basepath, dirnames{i}));
        [~, basename] = fileparts(filepath);
        
        if ~isempty(state)
            data = sr.states.fr{state};
            tstamps = sr.states.tstamps{state};
        else
            data = sr.strd;
            tstamps = sr.tstamps;
        end
        
        % find idx of CNO on/off
        onSession = 'lh70_201015_1815';
        offSession = 'lh70_201019_1831'; 
        if strcmp(basename, onSession)
            CNOidx(1) = length(msr);
        elseif strcmp(basename, offSession)
            CNOidx(2) = length(msr);
        else
            CNOidx = [0 0];
        end       
        
        % arrange time indices
        recStart = split(basename, '_');
        recStart = recStart{end};
        if numel(recStart) == 6
            tformat = 'HHmmss';
        elseif numel(recStart) == 4
            tformat = 'HHmm';
        end
        recStart = datetime(recStart, 'InputFormat', tformat);
        t = '120000';
        t = datetime(t, 'InputFormat', tformat);
        if t <= recStart
            t = t + hours(12);
        end
        s = seconds(t - recStart);
        [~, idx] = min(abs(s - tstamps));
        xidx(i) = idx + length(msr);      
        xname{i} = [char(sessionDate(i)) '_' datestr(t, 'HHMM')];
        xname{i} = datestr(t, 'HHMM');
        
        msr = [msr, data];
    end
    
    % smooth data
    gk = gausswin(600);
    gk = gk / sum(gk);
    for i = 1 : size(msr, 1)
        msr(i, :) = conv(msr(i, :), gk, 'same');
    end
    
    fh = figure;
    stdshade(msr(grp, :), '0.3', 'k')
    hold on
    yLimit = ylim;
    patch([CNOidx CNOidx(2) CNOidx(1)],...
        [yLimit(1) yLimit(1) yLimit(2) yLimit(2)],...
        'b', 'FaceAlpha', 0.1)
    xticks(xidx)
    xticklabels(xname)
    xtickangle(45)
    axis tight
    ylabel('mean SR [Hz]')
    xlabel('Time [hhmm]')
    title(txt)
    
%     fh = figure;
%     stdshade(movstd(msr', 199)', 0.3, 'k')
%     hold on
%     yLimit = ylim;
%     patch([CNOidx CNOidx(2) CNOidx(1)],...
%         [yLimit(1) yLimit(1) yLimit(2) yLimit(2)],...
%         'b', 'FaceAlpha', 0.1)
%     xticks(xidx)
%     xticklabels(xname)
%     xtickangle(45)
%     axis tight
%     ylabel('STD')
%     xlabel('Time [hhmm]')
end