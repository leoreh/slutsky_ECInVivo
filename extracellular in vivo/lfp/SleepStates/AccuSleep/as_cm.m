function [statePrecision, stateRecall] = as_cm(labels1, labels2, varargin)

% calculate and plot confusion matrix between labels1 (gold standard) and
% labels2. saves figure in graphics folder
%
% INPUT:
%   labels1         numeric. gold standard (e.g. manual scoring)
%   labels2         numeric. data for comparison (e.g. ntework output)
%   basename        string. name to assign figure
%   graphics        logical. plot confusion chart {true}
%   saveFig         logical. save figure {true}
%
% DEPENDENCIES:
%   as_loadConfig 
%   setMatlabGraphics
%
% 08 jun 21 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveFig', true, @islogical);

parse(p, varargin{:})
graphics        = p.Results.graphics;
saveFig         = p.Results.saveFig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = pwd;
[~, basename] = fileparts(basepath);

% load params from config file
[~, cfg_names, ~] = as_loadConfig([]);
nstates = length(cfg_names);

% select relavent labels 
idxLabels = labels2 < nstates & labels1 < nstates;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare functions
% precision (column): when predicts yes, how often is it correct? TP / (TP
% + FP). recall (row): when actually yes, how often does it predict yes? TP
% / (TP + FN)
precision = @(cm) diag(cm)./sum(cm, 2);
recall = @(confusionMat) diag(confusionMat) ./ sum(confusionMat, 1)';

% calc confusion matrix
cm = confusionmat(labels1(idxLabels), labels2(idxLabels));

statePrecision = precision(cm);
stateRecall = recall(cm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setMatlabGraphics(true)
if graphics
    fh = figure;
    fh.Position = [500 200 900 700];
    cmh = confusionchart(cm, cfg_names(1 : nstates - 1), 'ColumnSummary',...
        'column-normalized', 'RowSummary', 'row-normalized',...
        'title', 'State Classification Confusion Matrix', 'Normalization',...
        'total-normalized');
    sortClasses(cmh, cfg_names(1 : nstates - 1))
    
    if saveFig
        figpath = fullfile('graphics', 'sleepState');
        mkdir(figpath)
        figname = fullfile(figpath, sprintf('%s_CM', basename));
        export_fig(figname, '-tif', '-transparent', '-r300')
    end
end

end

% EOF


