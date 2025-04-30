function otl = get_otl(data, varargin)

% Generalized function to detect outliers accroding to distance in feature
% space. so far used for spectrogram and PSD data.
%
% INPUT
%   data        numeric mat of cases x variables (e.g., spectrogram or PSD bouts)
%   thrFactor   Threshold factor for outlier detection {5}.
%   thrMet      string. method to detect outliers {'mean'}
%   nfet        Number of principal components to use {3}.
%   nrep        Number of iterations {3}.
%   graphics    Logical. Plot detection results {false}.
%
% OUTPUT
%   otl         Struct containing outlier indices and detection info.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse inputs
p = inputParser;
addParameter(p, 'thrFactor', 1, @isnumeric);
addParameter(p, 'thrMet', 'mean', @ischar);
addParameter(p, 'nfet', 3, @isnumeric);
addParameter(p, 'nrep', 3, @isnumeric);
addParameter(p, 'graphics', false, @islogical);
parse(p, varargin{:});

thrFactor       = p.Results.thrFactor;
thrMet          = p.Results.thrMet;
nfet            = p.Results.nfet;
nrep            = p.Results.nrep;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterative Outlier Detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
tmp_data = data;
idx_otl = [];                    % outlier indices to original data mat (global)
idx_tmp = [1 : size(data, 1)]';  % for tracking original indices
% idx_iter                       % outlier indices to updated data after removal of outliers from previous iteration

% iterate
for irep = 1 : nrep

    % update thereshold
    thrVal = thrFactor * (1.5 ^ (irep - 1));
    % thrVal = thrFactor;

    % detect
    iter(irep) = find_otl(tmp_data, thrVal, nfet, thrMet);

    % break if no outliers detected
    if isempty(iter(irep).idx_iter)
        iter(irep) = [];
        break
    end

    % Track outliers globally
    iter(irep).idx_otl = idx_tmp(iter(irep).idx_iter);
    idx_otl = unique([idx_otl; iter(irep).idx_otl]);

    % Remove outliers from tmp_data and update index map for next iteration
    tmp_data(iter(irep).idx_iter, :) = [];
    idx_tmp(iter(irep).idx_iter) = [];

end

% Organize results
otl.idx = idx_otl;
otl.boolean = false(size(data, 1), 1);
otl.boolean(otl.idx) = true;
otl.bouts = binary2bouts('vec', otl.boolean, 'minDur', [], 'maxDur', [],...
    'interDur', [], 'exclude', false, 'flgPrnt', false);
otl.iter = iter;
otl.info.calcTime = datetime("now");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    fh = figure;
    colormap('parula');  % Heatmap color scheme

    nrep = length(otl.iter);
    retained_idx = (1 : size(data, 1))';

    % normalize data for heat map
    norm_data = zscore(data, 0, 1);  % Z-score normalization per column
    clipLimit = prctile(norm_data(:), [2, 98]);  % 2nd and 98th percentiles
    norm_data = min(max(norm_data, clipLimit(1)), clipLimit(2));

    for irep = 1 : nrep

        iter_data = otl.iter(irep);

        % Mean ± SEM with Outliers
        axh = subplot(nrep, 4, (irep - 1) * 4 + 1);
        retained_idx = setdiff(1 : size(data, 1), iter_data.idx_otl)';
        mean_ret = mean(data(retained_idx, :), 1);
        sem_ret = std(data(retained_idx, :), 0, 1) / sqrt(length(retained_idx));
        fill([1:size(data, 2), size(data, 2):-1:1], ...
            [mean_ret + sem_ret, fliplr(mean_ret - sem_ret)], ...
            [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        hold on
        ph = plot(data(iter_data.idx_iter, :)', 'Color', [1, 0, 0, 0.4], 'LineWidth', 0.5);
        ph(2) = plot(mean_ret, 'k', 'LineWidth', 2);
        xlabel('Features');
        title('Mean ± SEM');
        xlim([1, size(data, 2)]);
        set(axh, 'xscale', 'log', 'yscale', 'log')
        lgdTxt = sprintf('Outliers n=%d', length(iter_data.idx_iter));
        legend(ph, {lgdTxt, 'Mean'}, 'Location', 'southwest')

        % Heatmap - Sort data by distance
        subplot(nrep, 4, (irep - 1) * 4 + 2);
        [~, sort_idx] = sort(iter_data.dists, 'descend');
        imagesc(norm_data(sort_idx, :));
        colorbar;
        title(sprintf('Heatmap'));
        ylabel('Cases (sorted)');
        xlabel('Features');

        % Scree Plot
        subplot(nrep, 4, (irep - 1) * 4 + 3);
        plot(cumsum(iter_data.expl), 'LineWidth', 2);
        xlabel('PC');
        ylabel('Cum. Var. (%)');
        title('Scree Plot');
        hold on;
        xlim([0 10]);
        bh = bar(iter_data.expl);
        bh.FaceColor = 'flat';
        bh.CData(1 : nfet, :) = repmat([0 0 0], nfet, 1);
        bh.CData(nfet + 1 : end, :) = repmat([1 0 0], size(bh.CData, 1) - nfet, 1);
        ylabel('Explained Variance');
        xlabel('Principal Component');

        % PCA Projection
        axh = subplot(nrep, 4, (irep - 1) * 4 + 4);
        pcMat = iter_data.pc;
        badIdx = iter_data.idx_iter;
        % sh = scatter3(axh, pcMat(:, 1), pcMat(:, 2), pcMat(:, 3), ...
        %     20, 'k', 'filled');
        % sh.MarkerFaceAlpha = 0.4;
        hold on;
        scatter3(axh, pcMat(badIdx, 1), pcMat(badIdx, 2), pcMat(badIdx, 3), ...
            50, 'd', 'LineWidth', 1, 'MarkerEdgeColor', 'k');
        view(-25, 25);
        grid minor;
        xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
        title(axh, 'Feature Space');
    end
    sgtitle('Iterative Outlier Detection');
end


end

% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function iter_data = find_otl(tmp_data, thrVal, nfet, thrMet)

% detect outliers from data mat

% PCA
[~, pc, ~, ~, expl] = pca(tmp_data, 'NumComponents', nfet);

% Calculate distances
% dists = pc(:, 1);  % Use PC1 initially
dists = mahal(pc, pc);  % Final pass uses Mahalanobis distance

% Find outlier
switch thrMet
    case 'mad'
        thr = median(dists) + thrVal * mad(dists, 1);
        idx_iter = find(dists > thr);
    case 'mean'
        thr = mean(dists) + std(dists) * thrVal;
        idx_iter = find(dists > thr);
end

% Store iteration details
iter_data.pc = pc;
iter_data.expl = expl;
iter_data.dists = dists;
iter_data.thr = thr;
iter_data.idx_iter = idx_iter;
iter_data.idx_otl = [];

end
