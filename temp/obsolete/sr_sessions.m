% sr_sessions
% spike rate per tetrode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
forceL = false;
forceA = false;

% full path and name to xls file with session metadata
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
mname = 'lh91';     % mouse name

% column name in xls sheet where dirnames exist
colName = 'Session';

% string array of variables to load
vars = ["session.mat";...
    "SleepState.states.mat";...
    ".datInfo";
    ".sr.mat";...
    "AccuSleep_states"];

% column name of logical values for each session. only if true than session
% will be loaded. can be a string array and than all conditions must be
% met.
pcond = ["tempflag"];

% same but imposes a negative condition
ncond = ["fepsp"];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%
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
        [~, basename] = fileparts(basepath);
        cd(basepath)
        
        % params
        session = CE_sessionTemplate(pwd, 'viaGUI', false,...
            'force', true, 'saveVar', true);
        nchans = session.extracellular.nChannels;
        fs = session.extracellular.sr;
        spkgrp = session.extracellular.spikeGroups.channels;
        
        % pre-process dat
        %         mapch = 1 : 15;
        %         rmvch = [];
        %         datInfo = preprocDat('basepath', basepath, 'fname', [basename, '.dat'], 'mapch', mapch,...
        %             'rmvch', rmvch, 'nchans', nchans, 'saveVar', true,...
        %             'chunksize', 5e6, 'precision', 'int16', 'bkup', true,...
        %             'clip', round([406 * fs * 60 438 * fs * 60]), 'rmvArt', false);
        
        % create lfp
        %         LFPfromDat('basepath', basepath, 'cf', 450, 'chunksize', 5e6,...
        %             'nchans', nchans, 'fsOut', 1250,...
        %             'fsIn', fs)
        %
        % sleep states
        [EMG, EEG, sigInfo] = as_prepSig([basename, '.lfp'], [basename, '.emg.dat'],...
            'eegCh', 1, 'emgCh', 1, 'saveVar', true, 'emgNchans', 1,...
            'inspectSig', true, 'forceLoad', true, 'eegFs', 1250, 'emgFs', 6103.515625);
        %         AccuSleep_viewer(EEG, EMG, 1250, 1, [], [])
        
        % detect spikes
        %         [spktimes, ~] = spktimesWh('basepath', basepath, 'fs', fs, 'nchans', nchans,...
        %             'spkgrp', spkgrp, 'saveVar', true, 'saveWh', true,...
        %             'graphics', false, 'force', true);
        %
        %         % create ns files for sorting
        %         dur = [];
        %         t = [];
        %         spktimes2ns('basepath', basepath, 'fs', fs,...
        %             'nchans', nchans, 'spkgrp', spkgrp, 'mkClu', true,...
        %             'dur', dur, 't', t, 'psamp', [], 'grps', [1 : length(spkgrp)],...
        %             'spkFile', 'temp_wh');
        %
        %         % spike rate per tetrode. note that using firingRate requires
        %         % special care becasue spktimes is given in samples and not seconds
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

close all
tet = [1 : 4];          % which tetrodes to plot
state = [];             % [] - all; 1 - awake; 2 - NREM; 3 - REM
Y = [0 250];            % ylim
p1 = 1;                 % firing rate vs. time, one fig per session
p2 = 0;                 % mfr across sessions, one fig
p3 = 0;                 % spike rate vs. time, one fig for all sessions
plotStyle = 'bar';      % for p2. can by 'bar' or 'box'
clr = ['bbbbbkkkkkkkkrrrrrrrr'];        % color sessions
clr = ['bbbbbkkkkkkkkkkkkrrrrr'];        % color sessions
sessions = 1 : nsessions;
% sessions = 11;
[nsub] = numSubplots(length(sessions));

k = 1;
for i = sessions
    
    session = varArray{i, 1}.session;
    %     ss = varArray{i, 2}.SleepState;
    if ~isempty(varArray{i, 3})
        datInfo = varArray{i, 3}.datInfo;
    end
    sr = varArray{i, 4}.fr;
    if ~isempty(varArray{i, 5})
        ss = varArray{i, 5}.ss;
    else
        ss = [];
    end
    
    tstamps = sr.tstamps ;
    if isfield(datInfo, 'nsamps')
        nsamps = cumsum(datInfo.nsamps);
    elseif isfield(datInfo, 'blockduration')
        nsamps = cumsum(datInfo.blockduration) * fs;    % tdt legacy
    else
        nsamps = 120 * 60 * fs;
    end
    
    % mfr in selected state across sessions (mean +- std)
    %     states = {ss.ints.WAKEstate, ss.ints.NREMstate, ss.ints.REMstate};
    if ~isempty(ss)
        states = ss.boutTimes;
    else
        states = [];
    end
    if ~isempty(state)
        mfr{i} = mean(sr.states.fr{state}(tet, :), 2);
        txt = sprintf('MSR of T#%s in %s', num2str(tet), sr.states.statenames{state});
    else
        mfr{i} = mean(sr.strd(tet, :), 2);
        txt = sprintf('MSR of T#%s', num2str(tet));
    end
    
    % firing rate vs. time. 1 fig per session
    if p1
        if ~exist('fh1', 'var')
            fh1 = figure;
        end
        subplot(nsub(1), nsub(2), k)
        plot(sr.tstamps / 60, (sr.strd(tet, :))', 'LineWidth', 1.5)
        hold on
        % medata = median((data(grp, :)), 1);
        % plot(tstamps, medata, 'k', 'LineWidth', 5)
        % stdshade(sr.strd(grp, :), 0.3, 'k', tstamps / 60)
%         for ii = 1 : length(nsamps) - 1
%             plot([nsamps(ii) nsamps(ii)] / fs / 60, Y, '--k',...
%                 'LineWidth', 2)
%         end
        axis tight
        if ~isempty(states)
            fill([states{2} fliplr(states{2})]' / 60, [Y(1) Y(1) Y(2) Y(2)],...
                'b', 'FaceAlpha', 0.15,  'EdgeAlpha', 0);
        end
        ylim(Y)
        ylabel('Spike Rate [Hz]')
        xlabel('Time [m]')
        title(dirnames{i})
        if i == sessions(1)
            legend(split(num2str(tet)))
        end
        if saveFig
            figname = sprintf('%s_SRvsTime', sessionDate{i})
            figname = fullfile(mousepath, figname);
            % print(fh, figname, '-dpdf', '-bestfit', '-painters');
            export_fig(figname, '-tif', '-transparent', '-r300')
        end
    end
    k = k + 1;
end

if p2
    fh2 = figure;
    mfrmat = cell2nanmat(mfr);
    switch plotStyle
        case 'bar'
            plot([1 : nsessions], mfrmat, 'LineWidth', 1.5)
            hold on
            bh = bar(nanmean(mfrmat), 'FaceAlpha', 0.1, 'FaceColor', 'k');
            hold on
            % errorbar(1 : nsessions, nanmean(mfrmat), nanstd(mfrmat), 'k')
            idx = [find(clr == 'k', 1) find(clr == 'r', 1)] - 0.5;
            yLimit = ylim;
            patch([idx idx(2) idx(1)],...
                [yLimit(1) yLimit(1) yLimit(2) yLimit(2)],...
                'b', 'FaceAlpha', 0.1)
            legend(split(num2str(tet)))
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
        figname = fullfile(mousepath, ['sessionSummary of ' txt]);
        % print(fh, figname, '-dpdf', '-bestfit', '-painters');
        export_fig(figname, '-tif', '-transparent', '-r300')
    end
end

% spike rate vs. time across all sessions (1 fig)
% state = [];
onSession = 'lh86_210226_184500';
offSession = 'lh86_210228_190000';
patchIdx = [0 0];
p4 = 0;
if p4
    msr = [];
    tlast = 0;
    tss = [];
    NREMstates = [];
    prev_t = 0;
    for i = 1 : nsessions
        session = varArray{i, 1}.session;
        % ss = varArray{i, 2}.SleepState;
        if ~isempty(varArray{i, 3})
            datInfo = varArray{i, 3}.datInfo;
        end
        sr = varArray{i, 4}.fr;
        % st = varArray{i, 5}.spktimes;
        basepath = char(fullfile(mousepath, dirnames{i}));
        [~, basename] = fileparts(basepath);
        
        if ~isempty(state)
            data = sr.states.fr{state};
            tstamps = sr.states.tstamps{state};
        else
            data = sr.strd;
            tstamps = sr.tstamps;
            % tss = [tss; ss.ints.NREMstate + offset];
        end
        
        % find idx of CNO on/off
        if strcmp(basename, onSession)
            patchIdx(1) = length(msr);
        elseif strcmp(basename, offSession)
            patchIdx(2) = length(msr);
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
        t = '110000';
        t = datetime(t, 'InputFormat', tformat);
        if t <= recStart
            t = t + hours(12);
        end
        s = seconds(t - recStart);
        [~, idx] = min(abs(s - tstamps));
        xidx(i) = idx + length(msr);
        xname{i} = [char(sessionDate(i)) '_' datestr(t, 'HHMM')];
        
        msr = [msr, data];
        tsr{i} = tstamps + tlast;
        tlast = tsr{i}(end);
    end
    
    % smooth data
    smf = [1];
    if smf > 2
        gk = gausswin(smf);
        gk = gk / sum(gk);
        for i = 1 : size(msr, 1)
            msr(i, :) = conv(msr(i, :), gk, 'same');
        end
    end
    
    fh = figure;
    plot(msr(tet, :)', 'LineWidth', 1)
    hold on
    %     stdshade(msr(grp, :), '0.3', 'k')
    yLimit = ylim;
    patch([patchIdx patchIdx(2) patchIdx(1)],...
        [yLimit(1) yLimit(1) yLimit(2) yLimit(2)],...
        'r', 'FaceAlpha', 0.2)
    if ~isempty(state)
        fill([tss fliplr(tss)]' / 60, [yLimit(1) yLimit(1) yLimit(2) yLimit(2)],...
            'b', 'FaceAlpha', 0.15,  'EdgeAlpha', 0);
    end
    xticks(xidx)
    xticklabels(xname)
    xtickangle(45)
    axis tight
    ylabel('mean SR [Hz]')
    xlabel('Time [hhmm]')
    %     title(txt)
    if saveFig
        figname = fullfile(mousepath, txt);
        % print(fh, figname, '-dpdf', '-bestfit', '-painters');
        export_fig(figname, '-tif', '-transparent', '-r300')
    end
end