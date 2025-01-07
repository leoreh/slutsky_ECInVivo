function [tmp_idx, iter] = find_otl(tmp_data, thrVal)

        % detect outliers from data mat

        % Perform PCA
        [~, pc, ~, ~, expl] = pca(tmp_data, 'NumComponents', numPC);

        % Calculate distances
        dists = pc(:, 1);  % Use PC1 initially
        dists = mahal(pc, pc);  % Final pass uses Mahalanobis distance
        
        % Find outlier
        thr = median(dists) + thrVal * mad(dists, 1);
        tmp_idx = find(dists > thr);  % Outlier indices for current iteration

        if isempty(tmp_idx)
            return
        end

        % Store iteration details
        iter.pc = pc;
        iter.expl = expl;
        iter.dists = dists;
        iter.thr = thr;
        iter.outliers = tmp_idx;

    end