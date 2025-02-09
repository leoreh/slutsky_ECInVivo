function [totDur, prctDur, boutLen, timebins] = as_plotZT(varargin)

% receives a labels vector, divides it to states within timebins, and plots
% the bout destribution and duration. 
%
% INPUT:
%   basepath        char. path of recording
%   timebins        n x 2 numeric of time windows (indices to labels). 
%                   bouts of each state will be calculated for each time
%                   window separately. 
%   nbins           numeric. if timebins not specified, will divide labels
%                   to nbins timebins
%   sstates         numeric. idx of selected states. e.g. [1, 4, 5] will
%                   only plot wake, nrem and rem
%   ss              struct. see as_classify. if given then labels, minDur,
%                   and interDur will be extracted from ss
%   labels          numeric. see as_classify 
%   minDur          numeric. minimum duration of an bout. 
%                   if length(minDur) == nstates than a different minimum
%                   duration will be applied to each state
%   interDur        numeric. combine bouts separated by <= interDur
%   graphics        logical. 
% 
% OUTPUT
%
% DEPENDENCIES
% 
% TO DO LIST
%
% 12 jan 22 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'timebins', [], @isnumeric);
addParameter(p, 'nbins', [1], @isnumeric);
addParameter(p, 'sstates', [], @isnumeric);
addParameter(p, 'ss', []);
addParameter(p, 'labels', [], @isnumeric);
addParameter(p, 'minDur', [], @isnumeric);
addParameter(p, 'interDur', 4, @isnumeric);
addParameter(p, 'graphics', true, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
timebins        = p.Results.timebins;
nbins           = p.Results.nbins;
ss              = p.Results.ss;
sstates         = p.Results.sstates;
labels          = p.Results.labels;
minDur          = p.Results.minDur;
interDur        = p.Results.interDur;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file
[~, basename] = fileparts(basepath);
ssfile = fullfile(basepath, [basename, '.sleep_states.mat']);

% state params
cfg = as_loadConfig();
nstates = cfg.nstates;
if isempty(sstates)
    sstates = 1 : nstates;      % selected states (ignore bin)
end

% load data
if isempty(ss)
    if exist(ssfile, 'file')
        load(ssfile, 'ss')
    end
end
if isempty(labels)
    labels = ss.labels;
end

% arrange bout duration constraints
if isempty(minDur)
    minDur = ss.info.minDur;
end
if isempty(interDur)
    interDur = ss.info.interDur;
end
if length(minDur) == 1
    minDur = repmat(minDur, nstates, 1);
elseif length(minDur) ~= nstates
    error('minDur length is different than the number of states')
end
if length(interDur) == 1
    interDur = repmat(interDur, nstates, 1);
elseif length(interDur) ~= nstates
    error('interDur length is different than the number of states')
end

% create timebins
if isempty(timebins)
    timebins = n2chunks('n', length(labels), 'nchunks', nbins);
end
nbins = size(timebins, 1);
binLen = diff(timebins');

% calc bout stats according to timebins
for iwin = 1 : nbins
    labels_win = labels(timebins(iwin, 1) : timebins(iwin, 2));
    bouts(iwin) = as_bouts('labels', labels_win,...
        'minDur', minDur, 'interDur', interDur);
end

% organize bout stats
for iwin = 1 : nbins
    for istate = 1 : length(sstates)
        boutLen_temp{iwin, istate} = bouts(iwin).boutLen{sstates(istate)};
        if isempty(boutLen_temp{iwin, istate})
            boutLen_temp{iwin, istate} = 0;
        end
    end
    totDur(iwin, :) = bouts(iwin).totDur(sstates);
end
for istate = 1 : length(sstates)
    boutLen{istate} = cell2padmat(boutLen_temp(:, istate), 2);
end

% calc state duration [%]
for istate = 1 : length(sstates)
    prctDur(:, istate) = (totDur(:, istate) ./ binLen') * 100;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    setMatlabGraphics(true)
    alphaIdx = linspace(0.5, 1, nbins);
    fh = figure;
    
    for istate = 1 : length(sstates)
        
        % state duration
        subplot(2, length(sstates), istate)
        plot([1 : nbins], prctDur(:, istate),...
            'Color', cfg.colors{sstates(istate)}, 'LineWidth', 2)
        ylabel('State duration [min]')
        ylabel('State duration [%]')
        ylim([0 100])
        ax = gca;
        set(ax.YAxis, 'color', cfg.colors{sstates(istate)})
%         xticklabels(tbins_txt)
%         xtickangle(45)
        
        % bout length
        epMat = boutLen{istate};
        if ~all(epMat)
            continue
        end
        subplot(2, length(sstates), istate + length(sstates))
        boxplot(epMat, 'PlotStyle', 'traditional', 'Whisker', 100);
        bh = findobj(gca, 'Tag', 'Box');
        bh = flipud(bh);
        for ibox = 1 : length(bh)
            patch(get(bh(ibox), 'XData'), get(bh(ibox), 'YData'),...
                cfg.colors{sstates(istate)}, 'FaceAlpha', alphaIdx(ibox))
        end
        ylabel('Bout Length [log(s)]')
        set(gca, 'YScale', 'log')
        ylim([0 ceil(prctile(epMat(:), 99.99))])
        ax = gca;
        set(ax.YAxis, 'color', cfg.colors{sstates(istate)})
%         xticklabels(tbins_txt)
%         xtickangle(45)
        
    end
    sgtitle(basename)
    
    saveFig = true;
    if saveFig
        figpath = fullfile(basepath, 'graphics', 'sleepState');
        mkdir(figpath)
        figname = fullfile(figpath, [basename, '_ZTtimebins']);
        savefig(fh, figname, 'compact')
    end
end

end



% EOF