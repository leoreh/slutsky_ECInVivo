% fEPSP_Sessions

% rubust to missing sessions (replaced by nan) and allows for different
% stim intensities between sessions.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% should allow user to input varName or columnn index 
varName = 'fEPSP';
basepath = 'H:\Data\Processed\lh52';
sessionlist = 'sessionList.xlsx';       % must include extension

sessionInfo = readtable(fullfile(basepath, sessionlist));
col = strcmp(sessionInfo.Properties.VariableNames, varName);
dirnames = string(table2cell(sessionInfo(:, col)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : length(dirnames)
    filepath = fullfile(basepath, dirnames(i));
    if ~exist(filepath, 'dir')
        warning('%s does not exist, skipping...', filepath)
        continue
    end
    cd(filepath)
    filename = dir('*fepsp*');
    if length(filename) == 1
        f{i} = load(filename.name);
    else
        warning('no fepsp file in %s, skipping', filepath)
    end
end

% params
spkgrp = f{1}.fepsp.spkgrp;
ngrp = length(spkgrp);
nsessions = length(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rearrange data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intens = [];
for i = 1 : nsessions
    if isempty(f{i})
        continue
    end
    intens = sort(unique([intens, f{i}.fepsp.intens]));  
end

% 3d mat tetrodes x intensities x sessions
ampmat = nan(ngrp, length(intens), nsessions);

for i = 1 : nsessions
     if isempty(f{i})
        continue
     end
     samp = f{i}.fepsp.amp;
     sintens = f{i}.fepsp.intens;
     [~, ia] = intersect(intens, sintens);
     for ii = 1 : ngrp
        ampmat(ii, ia, i) = samp(ii, :);
     end
end

% amp during night (even) devided by values in day (odd).
% tetrodes x intensities x days
ndmat = ampmat(:, :, 2 : 2 : end) ./  ampmat(:, :, 1 : 2 : end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% one figure per intensity
for i = 1 : length(intens)
    figure
    plot(squeeze(ampmat(:, i, :))', 'LineWidth', 2)
    xlabel('Session')
    ylabel('Amplitude [mV]')
    title(sprintf('Stim Intensity %d uA', intens(i)))
    legend
end

% one figure per tetrode
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
    
% comparison night and day
si = 3;       % representative intensity
figure
plot(squeeze(ndmat(:, si, :))', 'LineWidth', 2)
hold on
plot([1 size(ndmat, 3)], [1 1], '--k')
axis tight
y = ylim;
ylim([0 y(2)])
xlabel('Time [days]')
ylabel('Ratio night / day')
legend(split(num2str(1 : ngrp)))
title('Ratio night / day')
