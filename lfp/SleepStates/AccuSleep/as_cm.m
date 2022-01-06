function [netPrecision, netRecall] = as_cm(labels1, labels2, varargin)

% calculate and plot confusion matrix between labels1 (gold standard) and
% labels2. if netScores are provided then will plot the effect of a
% confindance threshold on network performance. saves figure in graphics
% folder
%
% INPUT:
%   labels1         numeric. gold standard (e.g. manual scoring)
%   labels2         numeric. data for comparison (e.g. ntework output)
%   scores          numeric. network scores for each label.
%   graphics        logical. plot confusion chart {true}
%   saveFig         logical. save figure {true}
%
% DEPENDENCIES:
%   as_loadConfig
%   setMatlabGraphics
%
% 08 jun 21 LH      updates:
% 08 nov 21         score threshold 
% 02 jan 21         bug fix in score threshold

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'scores', [], @isnumeric);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveFig', true, @islogical);

parse(p, varargin{:})
scores          = p.Results.scores;
graphics        = p.Results.graphics;
saveFig         = p.Results.saveFig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = pwd;
[~, basename] = fileparts(basepath);

% load params from config file
cfg = as_loadConfig();
nstates = cfg.nstates;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare functions
% precision (column): when predicts yes, how often is it correct? TP / (TP
% + FP). recall (row): when actually yes, how often does it predict yes? TP
% / (TP + FN)
precision = @(cm) diag(cm)./sum(cm, 2);
recall = @(cm) diag(cm) ./ sum(cm, 1)';

% select relavent labels
idxLabels = labels2 <= nstates & labels1 <= nstates;
labels1 = labels1(idxLabels);
labels2 = labels2(idxLabels);
if ~isempty(scores)
    scores = scores(idxLabels, :);
end

% calc confusion matrix on entire data
cm = confusionmat(labels1, labels2);
netPrecision = precision(cm);
netRecall = recall(cm);

% simulate effect of thersholding the network score on performace
thr = 0 : 0.1 : 1;
thr_precision = zeros(length(thr), nstates);
thr_recall = zeros(length(thr), nstates);
for ithr = 1 : length(thr)
    scoreIdx = scores < thr(ithr);
    rmIdx = zeros(length(labels2), nstates);
    for istate = 1 : size(scores, 2)
        stateIdx = labels2 == istate;
        rmIdx(:, istate) = stateIdx & scoreIdx(:, istate);
        lostData(ithr, istate) = sum(rmIdx(:, istate)) / sum(stateIdx) * 100;
    end
    rmIdx = any(rmIdx, 2);
    labels1_thr = labels1(~rmIdx);
    labels2_thr = labels2(~rmIdx);
    cm_thr = confusionmat(labels1_thr, labels2_thr);
    temp_precision = precision(cm_thr);
    temp_recall = recall(cm_thr);
    
    % correct case where threholding removed all labels
    remainingStates = any(labels1_thr == [1 : 6]) | any(labels2_thr == [1 : 6]);
    
    thr_precision(ithr, remainingStates) = temp_precision;
    thr_recall(ithr, remainingStates) = temp_recall;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphics
    
    % ---------------------------------------------------------------------
    % figure of net performace
    setMatlabGraphics(true)
    fh = figure;
    fh.Position = [500 200 900 700];
    cmh = confusionchart(cm, cfg.names(1 : nstates), 'ColumnSummary',...
        'column-normalized', 'RowSummary', 'row-normalized',...
        'title', 'State Classification Confusion Matrix', 'Normalization',...
        'total-normalized');
    sortClasses(cmh, cfg.names(1 : nstates))
    
    if saveFig
        figpath = fullfile('graphics', 'sleepState');
        mkdir(figpath)
        figname = fullfile(figpath, sprintf('%s_cm', basename));
        export_fig(figname, '-tif', '-transparent', '-r300')
    end
        
    % ---------------------------------------------------------------------   
    % figure of net performace given threshold on net scores
    setMatlabGraphics(false)
    fh = figure;
    for istate = 1 : nstates
        
        subplot(2, size(scores, 2) / 2, istate)
        plot(thr, thr_precision(:, istate) * 100, 'k', 'LineWidth', 2)
        hold on
        plot(thr, thr_recall(:, istate) * 100, 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
        xlabel('Threshold')
        ylabel('Performance [%]')
        ylim([50 100])
        
        yyaxis right
        ph = plot(thr, lostData(:, istate), 'LineWidth', 2);
        ph.Color = cfg.colors{istate};
        ylabel('Data Lost [%]')
        set(gca, 'ycolor', cfg.colors{istate})
        ylim([0 100])
        
        set(gca, 'box', 'off', 'TickLength', [0 0])
        title(cfg.names{istate})
        if istate == 1
            legend({'Precision', 'Recall', 'DataLost'})
        end
    end
    
    if saveFig
        figpath = fullfile('graphics', 'sleepState');
        mkdir(figpath)
        figname = fullfile(figpath, sprintf('%s_cm_netScores', basename));
        export_fig(figname, '-tif', '-transparent', '-r300')
    end
end

end

% EOF


