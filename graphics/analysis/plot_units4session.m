
% plot a stacked bar plot of the number of units per session (top panel)
% and per spike group (bottom panel)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = 'lh99';
forceL = false;
if ~exist('varArray', 'var') || forceL
    [varArray, dirnames, mousepath] = getSessionVars('sortDir', false,...
        'dirnames', [], 'mname', mname);
end
nsessions = length(dirnames);
sessionIdx = 1 : nsessions;

% params
grp = [1 : 8];              % which tetrodes to plot
suFlag = 1;                 % plot only su or all units
% include only units with fr greater / lower than. 1st row RS 2nd row FS
frBoundries = [0.2 Inf; 0.2 Inf];    

[nsub] = numSubplots(length(sessionIdx));
setMatlabGraphics(false)
saveFig = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
maxUnits = 0;
units = cell(1, length(sessionIdx));

for isession = sessionIdx
    assignVars(varArray, isession)
    if isempty(varArray{isession, 3})
        units{isession}(:, 1) = zeros(length(grp), 1);
        units{isession}(:, 2) = zeros(length(grp), 1);
        continue
    end
    for igrp = unique(spikes.shankID)
        units{isession}(igrp, 1) = sum(selectUnits(spikes, cm, fr, suFlag, igrp, frBoundries, 'pyr'));
        units{isession}(igrp, 2) = sum(selectUnits(spikes, cm, fr, suFlag, igrp, frBoundries, 'int'));
    end
    maxUnits = max([maxUnits; sum(units{isession}, 2)]);
end

% arrange title names
for isession = 1 : nsessions
    sessionName{isession} = dirnames{isession}(length(mname) + 2 : end);
    basepath = char(fullfile(mousepath, dirnames{isession}));
    basepaths{isession} = fullfile(mousepath, dirnames{isession});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure;
k = 1;
for isession = sessionIdx
    subplot(nsub(1), nsub(2), k)
    bar(units{isession}, 'stacked')
    legend({"RS"; "FS"})
    xticks(1 : length(grp))
    title(sessionName{isession})
    xlabel('Spike Group')
    ylabel('No. Units')
    box off
    ylim([0 maxUnits])
    k = k + 1;
end

if saveFig
    mkdir(fullfile(basepath, 'graphics'))
    figname = fullfile(basepath, 'graphics', 'UnitsPerGrp');
    export_fig(figname, '-tif', '-transparent', '-r300')
end

% 
% clear units
% for isession = 1 : nsessions
%     assignVars(varArray, isession)
%     units(isession, 1) = sum(selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'pyr'));
%     units(isession, 2) = sum(selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'int'));
% end
% fh = figure;
% bar(units, 'stacked')
% legend({"RS"; "FS"})
% xticks(1 : nsessions)
% xticklabels(dirnames)
% xtickangle(45)
% title('Number of Units')
% xlabel('Session')
% ylabel('No. Units')
% box off

