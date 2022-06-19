function plot_hypnogram(varargin)

% plot simple hypnogram as colored lines

% INPUT
%   labels          1 x n integers of state labels. see as_classify 
%   stateEpochs     cell of n x 2 mats. see as_classify 
%   sstates         selected states to mark. see below
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
addOptional(p, 'axh', []);

parse(p,varargin{:})
labels          = p.Results.labels;
stateEpochs     = p.Results.stateEpochs;
sstates         = p.Results.sstates;
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

% remove small states
%     labels(labels == 6) = 5;
%     labels(labels == 3) = 4;

% re-calc state epochs from labels
if isempty(stateEpochs)
    minDur = [5, 2, 2, 5, 2, 2];
    interDur = 10;
    [stateEpochs, ~] = as_epochs('labels', labels,...
        'minDur', minDur, 'interDur', interDur);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

hold on
for istate = 1 : length(sstates)
    sepochs = stateEpochs{sstates(istate)};
    if ~isempty(sepochs)
        plot(axh, sepochs', 1 * ones(size(sepochs))',...
            'color', cfg.colors{sstates(istate)}, 'LineWidth', 25)
    end
end
set(gca, 'ytick', [])
set(gca, 'YColor', 'none')

end 

% EOF