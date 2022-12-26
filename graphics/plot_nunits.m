function plot_nunits(varargin)

% plot a stacked bar plot of the number of units per spike group
%
% INPUT:
%   basepaths       cell of chars to recording sessions
%   mname           char of mouse name. if basepaths is empty will get
%                   basepaths from the sessionList.xlsx file
%   saveFig         logical
%
% 26 12 22 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepaths', {}, @iscell);
addParameter(p, 'mname', '', @ischar);
addParameter(p, 'saveFig', true, @islogical);

parse(p, varargin{:})
basepaths       = p.Results.basepaths;
mname           = p.Results.mname;
saveFig         = p.Results.saveFig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data
varsFile = ["units"];
varsName = ["units"];
if isempty(basepaths)
    [v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', ["tempflag"]);
else
    [v, ~] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', ["tempflag"]);
end
nfiles = length(basepaths);
for ifile = 1 : nfiles
    [mousepath, basenames{ifile}] = fileparts(basepaths{ifile});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMatlabGraphics(false)
clr = [0.4 0.7 1; 0.95 0.4 0.4; 0.66 0.4 0.66];
fh = figure;

if nfiles == 1
    figpath = fullfile(basepaths{1}, 'graphics');
    figname = fullfile(figpath, [basenames{1}, '_nunits']);

    nunits = v(1).units.nunits;
    
    subplot(1, 2, 1)
    bh = bar(nunits', 'stacked');
    bh(1).FaceColor = clr(1, :);
    bh(2).FaceColor = clr(2, :);
    bh(3).FaceColor = clr(3, :);
    legend({'RS', 'FS', 'Other'})
    xlabel('Spike Group')
    ylabel('No. Units')
    box off

    subplot(1, 2, 2)
    bh = bar(sum(nunits', 1), 'FaceColor', 'flat');
    bh.CData = clr;
    ylabel('No. Units')
    xticklabels({'RS', 'FS', 'Other'})

    sgtitle(basenames{1})
    
else
    figpath = fullfile(mousepath, 'graphics');
    figname = fullfile(figpath, [mname, '_nunits']);
    
    % organize matrix [session (rows) x tetrode (columns)]
    units = catfields([v(:).units]);
    units = units.nunits;
    ngrps = size(units, 2);
    clear rs fs other
    for igrp = 1 : ngrps
        rs(:, igrp) = units([1 : 3 : end], igrp);
        fs(:, igrp) = units([2 : 3 : end], igrp);
        other(:, igrp) = units([3 : 3 : end], igrp);
    end
    
    subplot(1, 2, 1)
    plot(rs)
    hold on
    plot(sum(rs, 2), 'Color', 'k', 'LineWidth', 2)
    xticks([1 : ngrps])
    xticklabels(basenames)
    xlabel('session')
    ylabel('RS units')
    legend([split(num2str([1 : ngrps])); 'Total'])

    subplot(1, 2, 2)
    plot(fs)
    hold on
    plot(sum(fs, 2), 'Color', 'k', 'LineWidth', 2)
    xticks([1 : ngrps])
    xticklabels(basenames)
    xlabel('session')
    ylabel('FS units')
end

if saveFig
    mkdir(figpath)
    export_fig(figname, '-tif', '-transparent', '-r300')
end

end