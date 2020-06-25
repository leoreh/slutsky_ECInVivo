% fEPSP_Sessions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% should allow user to input varName or columnn index 
varName = 'fEPSP';
basepath = 'E:\Data\Processed\lh52';
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

spkgrp = f{1}.fepsp.spkgrp;
ngrp = length(spkgrp);
nsessions = length(f);
%%% how do you want to handle missing sessions? ignore or leave space?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rearrange data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create 3d mat of tetrodes x intensities x sessions. 
intens = [];
for i = 1 : nsessions
    if isempty(f{i})
        continue
    end
    intens = sort(unique([intens, f{i}.fepsp.intens]));  
end

amp = nan(ngrp, length(intens), nsessions);

for i = 1 : nsessions
     if isempty(f{i})
        continue
     end
     samp = f{i}.fepsp.amp;
     sintens = f{i}.fepsp.intens;
     [~, ia] = intersect(intens, sintens);
     for ii = 1 : ngrp
        amp(ii, ia, i) = samp(ii, :);
     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
% one figure per intensity
for i = 1 : length(intens)
    figure
    plot(squeeze(amp(:, i, :))')
    xlabel('Session')
    ylabel('Amplitude [mV]')
    title(sprintf('Stim Intensity %d uA', intens(i)))
    legend
end

% one figure per tetrode
for i = 1 : ngrp
    figure
    plot(squeeze(amp(i, :, :))')
    xlabel('Session')
    ylabel('Amplitude [mV]')
    title(sprintf('T%d', i))
    legend(split(num2str(intens)))
end
    

