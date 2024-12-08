function plot_hypnogram(varargin)

% plot simple hypnogram as colored lines 

% INPUT
%   labels          1 x n integers of state labels. see as_classify 
%   stateEpochs     cell of n x 2 mats. see as_classify 
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
addOptional(p, 'stateEpochs', [], @iscell);
addOptional(p, 'sstates', [], @isnumeric);
addOptional(p, 'clr', []);
addOptional(p, 'yshift', 1, @isnumeric);
addOptional(p, 'lWidth', 20, @isnumeric);
addOptional(p, 'axh', []);

parse(p,varargin{:})
labels          = p.Results.labels;
stateEpochs     = p.Results.stateEpochs;
sstates         = p.Results.sstates;
clr             = p.Results.clr;
yshift          = p.Results.yshift;
lWidth          = p.Results.lWidth;
axh             = p.Results.axh;

if isempty(stateEpochs) && isempty(labels)
    error('must input stateEpochs or labels')
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
% merging of nearby epochs. for example, if a certain bin belongs both to
% nrem and qw, setting state 4 (nrem) after state 2 (qw) will show the bin
% as nrem
if isempty(sstates)
    sstates = 1 : cfg.nstates;
end

% colors for states
if isempty(clr)
    clr = cfg.colors(sstates);
end

% re-calc state epochs from labels
if isempty(stateEpochs)
    [stateEpochs, ~] = as_epochs('labels', labels,...
        'minDur', minDur, 'interDur', interDur, 'graphics', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

hold on
yLimit = ylim;
for istate = 1 : length(sstates)
    sepochs = stateEpochs{sstates(istate)};
    if ~isempty(sepochs)
        plot(axh, sepochs', yLimit(2) * yshift * ones(size(sepochs))',...
            'color', clr{istate}, 'LineWidth', lWidth)
    end
end
set(gca, 'ytick', [])
set(gca, 'YColor', 'none')
yLimit = ylim;
ylim([yLimit(1), yLimit(2) * yshift])

end 

% EOF