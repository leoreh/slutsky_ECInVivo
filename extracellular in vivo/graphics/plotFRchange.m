function f = plotFRchange(mat, varargin)

% recieves a matrix and plots each row as a line. can also plot the mean +
% std across columns. labels for the x axis (columns) can be defined. for
% example, can be used to plot the change in firing rate under different
% conditions and/or time windows.

% INPUT
% required:
%   mat       matrix of struct (see calcFR)
% optional:
%   avg         plot mean {1} or not (0)
%   labels      array of strings. one for each column of avgfr
%   basepath    recording session path {pwd} to save figure
%   saveFig     save figure (1) or not {0}
% 
% 12 jan 18 LH.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'avg', true, @islogical);
addOptional(p, 'labels', [], @iscell);
addOptional(p, 'basepath', pwd);
addOptional(p, 'saveFig', false, @islogical);

parse(p,varargin{:})
avg = p.Results.avg;
labels = p.Results.labels;
basepath = p.Results.basepath;
saveFig = p.Results.saveFig;

[~, ncol] = size(mat);
if ncol < 2
    error('mat must contain more than one column')
end
if ~isempty(labels) && length(labels) ~= ncol
    error('the number of labels must be equal to the number of columns in mat')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure;
xpoints = 1 : ncol;
line(xpoints, mat, 'Color', [0.5 0.5 0.5])
hold on
line(xpoints, [mean(mat)], 'Color', 'k', 'LineWidth', 5)
axis tight
xlim([xpoints(1) - 0.2, xpoints(end) + 0.2])
xlabel('Treatment')
ylabel('Norm. Firing Rate')
% ylim([0 5])
if ~isempty(labels)
    ax = gca;
    ax.XTick = xpoints;
    ax.XTickLabel = labels;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveFig
    filename = 'changeFR';
    savePdf(filename, basepath, f)
end

end

% EOF