function [totDur, epLen] = as_plotZT(varargin)

% receives a labels vector, divides it to states within timebins, and plots
% the epoch destribution and duration. 
%
% INPUT:
%   basepath        char. path of recording
%   timebins        n x 2 numeric of time windows (indices to labels). 
%                   epochs of each state will be calculated for each time
%                   window separately. 
%   nwin            numeric. if timebins not specified, will divide labels
%                   to nwin timebins
%   sstates         numeric. idx of selected states. e.g. [1, 4, 5] will
%                   only plot wake, nrem and rem
%   ss              struct. see as_classify. if given then labels, minDur,
%                   and interDur will be extracted from ss
%   labels          numeric. see as_classify 
%   minDur          numeric. minimum duration of an epoch. 
%                   if length(minDur) == nstates than a different minimum
%                   duration will be applied to each state
%   interDur        numeric. combine epochs separated by <= interDur
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
addParameter(p, 'nwin', [1], @isnumeric);
addParameter(p, 'sstates', [], @isnumeric);
addParameter(p, 'ss', []);
addParameter(p, 'labels', [], @isnumeric);
addParameter(p, 'minDur', [10, 5, 5, 10, 5, 5], @isnumeric);
addParameter(p, 'interDur', 4, @isnumeric);
addParameter(p, 'graphics', true, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
timebins        = p.Results.timebins;
nwin            = p.Results.nwin;
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

if isempty(ss)
    if exist(ssfile, 'file')
        load(ssfile, 'ss')
    end
end
if isempty(labels)
    labels = ss.labels;
end
if isempty(minDur)
    minDur = ss.info.minDur;
end
if isempty(interDur)
    interDur = ss.info.interDur;
end


if length(minDur) == 1
    minDur = repamt(minDur, nstates, 1);
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
    timebins = n2chunks('n', length(labels), 'nchunks', nwin);
end
nwin = size(timebins, 1);

% calc epoch stats according to timebins
for iwin = 1 : nwin
    labels_win = labels(timebins(iwin, 1) : timebins(iwin, 2));
    [stateEpochs, epochStats(iwin)] = as_epochs('labels', labels_win,...
        'minDur', minDur, 'interDur', interDur);
end

% organize epoch stats
for iwin = 1 : nwin
    for istate = 1 : length(sstates)
        epLen_temp{iwin, istate} = epochStats(iwin).epLen{sstates(istate)};
    end
    totDur(iwin, :) = epochStats(iwin).totDur(sstates);
end
for istate = 1 : length(sstates)
    epLen{istate} = cell2nanmat(epLen_temp(:, istate));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    setMatlabGraphics(false)
    alphaIdx = linspace(0.5, 1, nwin);
    fh = figure;
    
    for istate = 1 : length(sstates)
        
        % state duration
        subplot(2, length(sstates), istate)
        epMat = epLen{istate};
        plot([1 : size(epMat, 2)], totDur(:, istate) / 60,...
            'Color', cfg.colors{sstates(istate)}, 'LineWidth', 2)
        ylabel('State duration [min]')
        ax = gca;
        set(ax.YAxis, 'color', cfg.colors{sstates(istate)})
%         xticklabels(tbins_txt)
%         xtickangle(45)
        
        % epoch length
        subplot(2, length(sstates), istate + length(sstates))
        boxplot(epMat, 'PlotStyle', 'traditional', 'Whisker', 100);
        bh = findobj(gca, 'Tag', 'Box');
        bh = flipud(bh);
        for ibox = 1 : length(bh)
            patch(get(bh(ibox), 'XData'), get(bh(ibox), 'YData'),...
                cfg.colors{sstates(istate)}, 'FaceAlpha', alphaIdx(ibox))
        end
        ylabel('Epoch Length [log(s)]')
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
        export_fig(figname, '-jpg', '-transparent', '-r300')
    end
end

end



% EOF