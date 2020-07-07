% fEPSP_Sessions

% organizes and plots fepsp from multiple sessions. gets structs and
% transforms them to matrices of vars (e.g. amp and wv) vs. time
% (sessions).rubust to missing sessions (replaced by nan) and allows for
% different stim intensities between sessions. compensates if arrays are
% not sorted (though fEPSPfromOE should sort by intensity)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
forceL = false;
forceA = false;

% should allow user to input varName or columnn index
varName = 'fEPSP';
basepath = 'H:\Data\Processed\lh52';
sessionlist = 'sessionList.xlsx';       % must include extension
fs = 20000;                             % can also be loaded from datInfo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if forceL
    
    % get files
    sessionInfo = readtable(fullfile(basepath, sessionlist));
    col = strcmp(sessionInfo.Properties.VariableNames, varName);
    dirnames = string(table2cell(sessionInfo(:, col)));
    
    for i = 1 : length(dirnames)
        filepath = fullfile(basepath, dirnames(i));
        if ~exist(filepath, 'dir')
            warning('%s does not exist, skipping...', filepath)
            continue
        end
        cd(filepath)
        
        filename = dir('*fepsp*');
        if length(filename) == 1
            load(filename.name);
            f{i} = fepsp;
        else
            warning('no fepsp file in %s, skipping', filepath)
        end
        if forceA
            intens = [100 50 150 200 250 300];
            fepsp = getfEPSPfromOE('basepath', char(filepath), 'fname', '', 'nchans', 27,...
                'spkgrp', f{i}.spkgrp, 'intens', intens, 'concat', false, 'saveVar', true,...
                'force', true, 'vis', 'off');
            f{i} = fepsp;
        end
    end
end

% params
spkgrp = f{1}.spkgrp;
ngrp = length(spkgrp);
nsessions = length(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rearrange data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% timestamps
tstamps = f{1}.t / fs * 1000;

% intensities
intens = [];
for i = 1 : nsessions
    if isempty(f{i})
        continue
    end
    intens = sort(unique([intens, f{i}.intens]));
end

% ampmat: 3d mat of amplitudes; tetrode x intensity x session
ampmat = nan(ngrp, length(intens), nsessions);
% wvmat: 3d mat of average waveforms; tetrode x session x sample (for 1
% selected intensity)
wvmat = nan(ngrp, nsessions, size(f{1}.wvavg, 3));
si = 250;
for i = 1 : nsessions
    if isempty(f{i})
        continue
    end
    samp = f{i}.amp;
    swv = f{i}.wvavg;
    sintens = f{i}.intens;
    [~, ia] = intersect(intens, sintens);
    [~, ib] = intersect(sintens, si);
    for ii = 1 : ngrp
        ampmat(ii, ia, i) = samp(ii, :);
        wvmat(ii, i, :) = swv(ii, ib, :);
    end
end

% samp: 2d mat for selected tetrode and intensity; rep x session
sg = 1;         % selected group
si = 250;       % selected intensity
clear samp
for i = 1 : nsessions
    if isempty(f{i})
        continue
    end
    sintens = f{i}.intens;
    [~, ib] = intersect(sintens, si);
    samp(i) = f{i}.ampcell(sg, ib);
end
samp = cell2nanmat(samp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% one figure per intensity
p = 1;
if p
    for i = 1 : length(intens)
        figure
        plot(squeeze(ampmat(:, i, :))', 'LineWidth', 2)
        xlabel('Session')
        ylabel('Amplitude [mV]')
        title(sprintf('Stim Intensity %d uA', intens(i)))
        legend
    end
end

% one figure per tetrode
p = 1;
if p
    for i = 1 : ngrp
        figure
        plot(squeeze(ampmat(i, :, :))', 'LineWidth', 2)
        axis tight
        y = ylim;
        ylim([0 y(2)])
        xticks(1 : 2 : nsessions)
        xticklabels(split(num2str(1 : nsessions / 2)))
        xlabel('Time [days]')
        ylabel('Amplitude [mV]')
        title(sprintf('T%d', i))
        legend(split(num2str(intens)))
        box off
    end
end

% waveform across sessions
p = 1;
ss = [4, 8, 16, 24];    % selected sessions
sg = [1, 7];
if p
    for i = sg
        figure
        plot(tstamps, squeeze(wvmat(i, ss, :))')
        axis tight
        xlabel('Time [ms]')
        ylabel('Amplitude [mV]')
        title(sprintf('T%d', i))
        box off
        legend(split(dirnames(ss)), 'Interpreter', 'none');
    end
end

% waveform across time within session
p = 0;
sg = 7;         % selected group
si = 250;       % selected intensity
ss = 1;         % selected session
if p
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
p = 1;
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

ndmat = ampmat(sg, ib, 2 : 2 : end) ./  ampmat(sg, ib, 1 : 2 : end);
if p
    figure
    [~, ib] = intersect(intens, si);
    plot(squeeze(mean(ndmat(:, :, :), 2))', 'LineWidth', 2)
    hold on
    plot([1 size(ndmat, 3)], [1 1], '--k')
    axis tight
    y = ylim;
    ylim([0 y(2)])
    xlabel('Time [days]')
    ylabel('Ratio night / day')
    legend(split(num2str(sg)))
    title(sprintf('night / day @%d uA', intens(ib)))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to prism
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
squeeze(ampmat(7, :, :));
squeeze(wvmat(7, :, :));

squeeze(mean(mean(ndmat(:, :, 1 : 3), 2), 1))
squeeze(mean(mean(ndmat(:, :, 5 : 10), 2), 1))
squeeze(mean(mean(ndmat(:, :, 12 : 13), 2), 1))

plot(squeeze(mean(ampmat([1, 4 : 8], [2 : 6], :), 1))')