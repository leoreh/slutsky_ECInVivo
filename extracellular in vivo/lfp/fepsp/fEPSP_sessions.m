% fEPSP_sessions

% organizes and plots fepsp from multiple sessions. gets structs and
% transforms them to matrices of vars (e.g. amp and wv) vs. time
% (sessions). rubust to missing sessions (replaced by nan) and allows for
% different stim intensities between sessions. compensates if arrays are
% not sorted (though fEPSPfromOE should sort by intensity)

forceA = false;
forceL = false;
saveFig = false;

% full path and name to xls file with session metadata
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
mname = 'lh76';     % mouse name

% column name in xls sheet where dirnames exist
colName = 'Session';

% string array of variables to load
vars = ["session.mat";...
    "fepsp"];

% column name of logical values for each session. only if true than session
% will be loaded. can be a string array and than all conditions must be
% met.
pcond = ["fepsp"; "tempflag"];

% same but imposes a negative condition
ncond = ["spktimes"];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('varArray', 'var') && ~forceL
    [varArray, dirnames, mousepath] = getSessionVars('vars', vars,...
        'pcond', pcond, 'ncond', ncond, 'sortDir', false, 'dirnames', [],...
        'xlsname', xlsname, 'mname', mname);
end
nsessions = length(dirnames);

% session info
basepath = char(fullfile(mousepath, dirnames(1)));
cd(basepath)
session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'force', true, 'saveVar', true);
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;
ngrp = length(spkgrp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if forceA
    close all
    for i = 1 : nsessions
        
        % file
        basepath = char(fullfile(mousepath, dirnames(i)));
        cd(basepath)
        
        fepsp = varArray{i, 2}.fepsp;
        
        % fepsp
        intens = [];
        fepsp = fEPSPfromDat('basepath', basepath, 'fname', '', 'nchans', nchans,...
            'spkgrp', spkgrp, 'intens', intens, 'saveVar', false,...
            'force', true, 'extension', 'dat', 'recSystem', 'oe',...
            'protocol', 'stp', 'anaflag', true, 'inspect', false, 'fsIn', fs,...
            'cf', 0);
        
%         fepsp = fEPSP_analysis('fepsp', fepsp, 'basepath', basepath,...
%             'force', true);
  
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rearrange data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% intensities throughout sessions
intens = [];
for i = 1 : nsessions
    if isempty(varArray{i, 2})
        continue
    end
    fepsp = varArray{i, 2}.fepsp;
    intens = sort(unique([intens, fepsp.intens]));
end

si = [];            % selected intensity [uA]. for stp if empty than maximum facilitation will bet taken
grp = 1;            % selected tetrode
dataVar = 'ampcell';
if exist('fepsp', 'var')
    protocol = fepsp.info.protocol;
else
    protocol = 'io';
end
switch protocol
    case 'io'
        % ampmat    3d mat of amplitudes; tetrode x intensity x session. data will
        %           be extrapolated to intensities not recorded in session.
        % iomat     2d mat of io curves; intensity x session
        % wvmat     3d mat of average waveforms; tetrode x session x sample
        %           for 1 selected intensity.
        % ampcell   array of amplitudes for selected intensity. each cell
        %           contains amps of all traces. if itensity not recorded will be
        %           extrapolated.
        ampmat = nan(ngrp, length(intens), nsessions);
        wvmat = nan(ngrp, nsessions, size(fepsp.traceAvg, 3));
        ampcell = cell(1, nsessions);
        iomat = nan(length(intens), nsessions);
        for i = 1 : nsessions
            fepsp = varArray{i, 2}.fepsp;
            sintens = sort(fepsp.intens);
            [~, ia] = intersect(sintens, si);
            [~, ib] = intersect(intens, si);
            [~, ic] = intersect(intens, sintens);
            ampmat(:, :, i) = [interp1(sintens, fepsp.amp(:, :)', intens, 'linear', 'extrap')]';
            if ~isempty(ia)
                for ii = 1 : ngrp
                    wvmat(ii, i, :) = fepsp.traceAvg(ii, ia, :);
                    ampcell{i} = fepsp.(dataVar){grp, ia};
                end
            else
                ampcell{i} = ampmat(grp, ib, i);
            end
            iomat(ic, i) = mean(cell2nanmat(fepsp.(dataVar)), 'omitnan');
        end
        
    case 'stp'
        wvmat = nan(ngrp, nsessions, size(fepsp.traceAvg, 3));
        ampmat = nan(nsessions, length(fepsp.ampNorm));
        for i = 1 : nsessions
            fepsp = varArray{i, 2}.fepsp;
            fac{i} = max(squeeze(fepsp.ampNorm(grp, : ,:))');
            if isempty(si)
                    [~, ia] = max(fac{i});
            end
            wvmat(:, i, :) = fepsp.traceAvg(:, ia, :);
            ampmat(i, :) = squeeze(fepsp.ampNorm(grp, ia, :));
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all

pathPieces = regexp(dirnames(:), '_', 'split'); % assumes filename structure: animal_date_time
sessionDate = [pathPieces{:}];
sessionDate = sessionDate(2 : 3 : end);
[nsub] = numSubplots(nsessions);

% waveform and box plot of amplitudes across sessions for selected
% intensity and tetrode. can select specific sessions for waveform
% ss = [1, 2, 3, 4, 8, 12];
ss = [1 : nsessions];
p = 1;
clr = ['kgggrrryyy'];        % must be sorted (i.e. g before k before r)
if p
    fh = figure;
    subplot(2, 1, 1)
    tstamps = [1 : size(wvmat, 3)] / fs * 1000;
    ph = plot(tstamps, squeeze(wvmat(grp, ss, :)), 'LineWidth', 2);
    clrRep = histc(clr, unique(clr));
    clear alphaIdx
    for i = 1 : length(clrRep)
        alphaIdx{i} = linspace(1, 0.3, clrRep(i));
    end
    alphaIdx = [alphaIdx{:}];
    if length(clr) == length(ph)
        for i = 1 : length(ph)
            ph(i).Color = clr(ss(i));
            ph(i).Color(4) = alphaIdx(ss(i));
        end
    end
    xlabel('Time [ms]')
    ylabel('Voltage [V]')
    legend(sessionDate(ss))
    box off
    
    subplot(2, 1, 2)
    switch protocol
        case 'io'
            ampmat = cell2nanmat(ampcell);
%                         boxplot(ampmat, 'PlotStyle', 'traditional');
            bar(nanmean(ampmat))
            bh = findobj(gca, 'Tag', 'Box');
            if length(bh) == length(clr)
                clr = fliplr(clr);
                alphaIdx = fliplr(alphaIdx);
                for i = 1 : length(bh)
                    patch(get(bh(i), 'XData'), get(bh(i), 'YData'),...
                        clr(i), 'FaceAlpha', alphaIdx(i))
                end
            end
            ylabel(dataVar)
            xticks(1 : nsessions)
            xticklabels(sessionDate)
            xtickangle(45)
            xlabel('Session')
        case 'stp'
            plot(ampmat')
            xticks(1 : length(fepsp.ampNorm))
            xlabel('Stimulus No.')
            ylabel('Normalized Amplitude')
    end
    y = ylim;
    %     ylim([0 y(2)]);
    
    box off
    suptitle([mname ' T#' num2str(grp) ' @ ' num2str(si) 'uA'])
    
    if saveFig
        figname = fullfile(mousepath, 'fepspSessions');
        % print(fh, figname, '-dpdf', '-bestfit', '-painters');
        export_fig(figname, '-tif', '-transparent', '-r300')
    end
end

% one figure per intensity across sessions
p = 0;
if p
    for i = 1 : length(intens)
        figure
        plot(squeeze(ampmat(:, i, :))', 'LineWidth', 2)
        xlabel('Session')
        ylabel('Amplitude [mV]')
        title(sprintf('Stim Intensity %d uA', intens(i)))
        legend(strsplit(num2str(1 : ngrp)))
        xticks(1 : nsessions)
        xticklabels(sessionDate)
        xtickangle(45)
        box off
    end
end

% io across sessions, one figure per selected grp
p = 0;
grp = 1;
if p
    for i = grp
        fh = figure;
        plot(iomat)
        axis tight
        xticklabels(split(num2str(intens)))
        xlabel('Intensity [uA]')
        ylabel('Amplitude [mV]')
        title(sprintf('T%d', i))
        legend(sessionDate)
        box off
    end
end

% one figure per tetrode
p = 0;
if p
    for i = 1 : ngrp
        figure
        plot(squeeze(ampmat(i, :, :))', 'LineWidth', 2)
        axis tight
        y = ylim;
        ylim([0 y(2)])
        xticks(1 : nsessions)
        xlabel('Session')
        ylabel('Amplitude [mV]')
        title(sprintf('T%d', i))
        legend(split(num2str(intens)))
        box off
    end
end

% waveform across time within session
p = 0;
if p
    sg = 7;         % selected group
    si = 250;       % selected intensity
    ss = 1;         % selected session
    for i = sg
        figure
        suptitle(sprintf('T%d @ %s', i, dirnames(ss)))
        sintens = f{ss}.intens;
        [~, ib] = intersect(sintens, si);
        swv = squeeze(mean(f{ss}.wv{i, ib}, 1));
        samp = f{ss}.ampcell{i, ib};
        
        subplot(1 ,2 ,1)
        plot(tstamps, swv)
        axis tight
        xlabel('Time [ms]')
        ylabel('Amplitude [mV]')
        legend
        box off
        
        subplot(1, 2, 2)
        plot(1 : length(samp), samp)
        xlabel('Stim #')
        ylabel('Amplitude [mV]')
        y = ylim;
        ylim([0 y(2)])
        box off
    end
end

% comparison night and day
% amp during night (even) devided by values in day (odd).
% tetrodes x intensities x days.
p = 0;
sg = [1, 4 : 8];    % excluded tetrodes not in ca1.
si = [200, 250, 300];    % reliable intensities
[~, ib] = intersect(intens, si);

% night = ampmat(sg, ib, 2 : 2 : end);
% night = mean(night(:, :, [1 : 3]), 3);
% day = ampmat(sg, ib, 1 : 2 : end);
% day = mean(day(:, :, [1 : 3]), 3);
% night ./ day
% night = ampmat(sg, ib, 2 : 2 : end);
% night = mean(night(:, :, [5 : 10]), 3);
% day = ampmat(sg, ib, 1 : 2 : end);
% day = mean(day(:, :, [5 : 10]), 3);
% night ./ day
% night = ampmat(sg, ib, 2 : 2 : end);
% night = mean(night(:, :, [12 : 13]), 3);
% day = ampmat(sg, ib, 1 : 2 : end);
% day = mean(day(:, :, [12 : 13]), 3);
% night ./ day

% ndmat = ampmat(sg, ib, 2 : 2 : end) ./  ampmat(sg, ib, 1 : 2 : end);
% if p
%     figure
%     [~, ib] = intersect(intens, si);
%     plot(squeeze(mean(ndmat(:, :, :), 2))', 'LineWidth', 2)
%     hold on
%     plot([1 size(ndmat, 3)], [1 1], '--k')
%     axis tight
%     y = ylim;
%     ylim([0 y(2)])
%     xlabel('Time [days]')
%     ylabel('Ratio night / day')
%     legend(split(num2str(sg)))
%     title(sprintf('night / day @%d uA', intens(ib)))
% end

% find range of stim artifact

fh = figure;
for i = 2 : nsessions
    
    % file
    basepath = char(fullfile(mousepath, dirnames(i)));
    cd(basepath)
    
    fepsp = varArray{i, 2}.fepsp;
    
    [~, stimidx] = min(abs(fepsp.tstamps - 0));
    margin = fepsp.info.fs * 1 / 1000;
    stimidx = [stimidx - margin : stimidx + margin];
    alphaIdx = linspace(1, 0.3, length(fepsp.intens));
    cmap = colormap(winter(length(fepsp.intens)));
    subplot(nsub(1), nsub(2), i)
    hold on
    clear stimArt
    for j = 1 : length(fepsp.intens)
        stimArt{j} = range(fepsp.traces{j}(stimidx, :));
        ph = plot(fepsp.ampcell{j}, stimArt{j}, 'LineStyle', 'none',...
            'Marker', 'o', 'MarkerFaceColor', cmap(j, :), 'MarkerSize', 8,...
            'MarkerEdgeColor', 'none');
        %         ph.Color = cmap(j)
    end
    xlabel('fEPSP Amplitude [mV]')
    ylabel('Sti. Artifact Amplitude [mV]')
    ylim([0 0.08])
    xlim([0. 0.6])
    
    [r, p] = corrcoef([fepsp.ampcell{:}], [stimArt{:}]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to prism
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% squeeze(ampmat(7, :, :));
% squeeze(wvmat(7, :, :));
%
% squeeze(mean(mean(ndmat(:, :, 1 : 3), 2), 1))
% squeeze(mean(mean(ndmat(:, :, 5 : 10), 2), 1))
% squeeze(mean(mean(ndmat(:, :, 12 : 13), 2), 1))
%
% plot(squeeze(mean(ampmat([1, 4 : 8], [2 : 6], :), 1))')