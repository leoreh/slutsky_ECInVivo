function [statePrecision, stateRecall] = as_cm(labels1, labels2, varargin)

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
[cfg_colors, cfg_names, ~] = as_loadConfig([]);
nstates = length(cfg_names);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare functions
% precision (column): when predicts yes, how often is it correct? TP / (TP
% + FP). recall (row): when actually yes, how often does it predict yes? TP
% / (TP + FN)
precision = @(cm) diag(cm)./sum(cm, 2);
recall = @(confusionMat) diag(confusionMat) ./ sum(confusionMat, 1)';

% select relavent labels
idxLabels = labels2 < nstates & labels1 < nstates;
labels1 = labels1(idxLabels);
labels2 = labels2(idxLabels);
if ~isempty(scores)
    scores = scores(idxLabels, :);
end

% calc confusion matrix
cm = confusionmat(labels1, labels2);
statePrecision = precision(cm);
stateRecall = recall(cm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphics
    
    if isempty(scores)
        setMatlabGraphics(true)
        fh = figure;
        fh.Position = [500 200 900 700];
        cmh = confusionchart(cm, cfg_names(1 : nstates - 1), 'ColumnSummary',...
            'column-normalized', 'RowSummary', 'row-normalized',...
            'title', 'State Classification Confusion Matrix', 'Normalization',...
            'total-normalized');
        sortClasses(cmh, cfg_names(1 : nstates - 1))
    
    else
        
        % simulate effect of score threshold on performace
        newLabels = labels2;
        thr = 0 : 0.1 : 1;       
        for ithr = 1 : length(thr)
            for istate = 1 : size(scores, 2)
                
                stateIdx = labels2 == istate;
                scoreIdx = scores(:, istate) < thr(ithr) & stateIdx;
                
                lostData(ithr, istate) = sum(scoreIdx) / sum(stateIdx) * 100;
                newLabels(scoreIdx) = nstates + 1;
            end
            
            % poor fix for when removes an entire state 
            tempCm = confusionmat(labels1, newLabels);
            tempPrecision = precision(cm);
            tempRecall = recall(cm);
%             [tempPrecision, tempRecall] =...
%                 as_cm(labels1, newLabels, 'graphics', false);
            if length(tempPrecision) == nstates - 1
                netPrecision(ithr, :) = tempPrecision;
                netRecall(ithr, :) = tempRecall;
            else
                netPrecision(ithr, :) = zeros(1, nstates - 1);
                netRecall(ithr, :) = zeros(1, nstates - 1);
            end
        end
        
        % plot performace vs. threshoold
        setMatlabGraphics(false)
        fh = figure;
        for istate = 1 : size(scores, 2)
            
            subplot(2, size(scores, 2) / 2, istate)
            plot(thr, netPrecision(:, istate) * 100, 'k', 'LineWidth', 2)
            hold on
            plot(thr, netRecall(:, istate) * 100, 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
            xlabel('Threshold')
            ylabel('Performance [%]')
            ylim([50 100])
            
            yyaxis right
            ph = plot(thr, lostData(:, istate), 'LineWidth', 2);
            ph.Color = cfg_colors{istate};
            ylabel('Data Lost [%]')
            set(gca, 'ycolor', cfg_colors{istate})
            ylim([0 100])
            
            set(gca, 'box', 'off', 'TickLength', [0 0])
            title(cfg_names{istate})
            if istate == 1
                legend({'Precision', 'Recall', 'DataLost'})
            end
        end
    end
    
    if saveFig
        figpath = fullfile('graphics', 'sleepState');
        mkdir(figpath)
        figname = fullfile(figpath, sprintf('%s_CM', basename));
        export_fig(figname, '-tif', '-transparent', '-r300')
    end
end

end

% EOF


