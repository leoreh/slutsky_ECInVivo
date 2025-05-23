function plot_hypnogram(varargin)

% plot simple hypnogram as colored lines 

% INPUT
%   labels          1 x n integers of state labels. see as_classify 
%   boutTimes       cell of n x 2 mats. see as_classify 
%   sstates         selected states to mark. see below
%   clr             cell of colors for each state. if empty will load cfg.colors
%   yshift          scalar. shift location of lines {1}
%   lWidth          scalar. width of lines to plot
%   axh             axis handle for plot
%
%   18 jun 22 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'labels', [], @isnumeric);
addOptional(p, 'boutTimes', [], @iscell);
addOptional(p, 'sstates', [], @isnumeric);
addOptional(p, 'clr', []);
addOptional(p, 'yshift', 1, @isnumeric);
addOptional(p, 'lWidth', 20, @isnumeric);
addOptional(p, 'axh', []);

parse(p,varargin{:})
labels          = p.Results.labels;
boutTimes       = p.Results.boutTimes;
sstates         = p.Results.sstates;
clr             = p.Results.clr;
yshift          = p.Results.yshift;
lWidth          = p.Results.lWidth;
axh             = p.Results.axh;

if isempty(boutTimes) && isempty(labels)
    error('must input boutTimes or labels')
end

if isempty(axh)
    fh = figure;
    axh = subplot(1, 1, 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

% state params
cfg = as_loadConfig();

% selected states. the order of sstates determines which state will be
% shown in case of overlap in sleepStates. This overlap occurs due to
% merging of nearby bouts. for example, if a certain bin belongs both to
% nrem and qw, setting state 4 (nrem) after state 2 (qw) will show the bin
% as nrem
if isempty(sstates)
    if isempty(labels)
        sstates = 1 : cfg.nstates;
    else
        sstates = unique(labels(~isnan(labels)));
    end
end

% colors for states
if isempty(clr)
    clr = cfg.colors(sstates);
end

% re-calc state bouts from labels
if isempty(boutTimes)
    bouts = as_bouts('labels', labels,...
        'minDur', 5, 'interDur', 3, 'graphics', false);
    boutTimes = bouts.times;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

hold on
yLimit = ylim;
for istate = 1 : length(sstates)
    sbouts = boutTimes{sstates(istate)};
    if ~isempty(sbouts)
        plot(axh, sbouts', yLimit(2) * yshift * ones(size(sbouts))',...
            'color', clr{istate}, 'LineWidth', lWidth)
    end
end
set(gca, 'ytick', [])
set(gca, 'YColor', 'none')
yLimit = ylim;
ylim([yLimit(1), yLimit(2) * yshift])

end 

% EOF