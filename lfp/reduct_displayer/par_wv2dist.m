function [dist_mat] = par_wv2dist(wv,dist_fun, assume_symmetry)
arguments
    wv (:,:) double {mustBeNumeric}
    dist_fun (1,1) function_handle = @(x,y) 1 - max(abs(xcorr(x,y,"normalized")));
    assume_symmetry (1,1) logical {mustBeNumericOrLogical} = true;
end

nWv = size(wv,1);
if assume_symmetry
    
    % collect linear index for loop to run on
    linIdxs = reshape(1:(nWv*nWv), nWv, nWv); % matrix of linear indices
    linIdxs = tril(linIdxs, -1);
    linIdxs(linIdxs == 0) = []; % Discard the indices that were set to 0 by TRIL
    
    % prealocate dist_mat
    dist_mat = nan(size(linIdxs));

    % use threads if it no pool exist, mainly for opening benefit
    if isempty(gcp("nocreate"))
        parpool("Threads",[1 floor(numel(linIdxs)*0.1)]);
    end

    % calculate distances
    parfor idx = 1:numel(linIdxs)
        % Convert back from a linear index to [i,j] co-ordinate pair
        [i,j] = ind2sub([nWv, nWv], linIdxs(idx));
        dist_mat(idx) = dist_fun(wv(i,:),wv(j,:));
    end
    dist_mat = squareform(dist_mat);
else
    % collect linear index for loop to run on
    linIdxs = reshape(1:(nWv*nWv), nWv, nWv); % matrix of linear indices
    linIdxs = linIdxs(:);

    % prealocate dist_mat
    dist_mat = nan(size(linIdxs));

    % use threads if it no pool exist, mainly for opening benefit
    if isempty(gcp("nocreate"))
        parpool("Threads",[1 floor(numel(linIdxs)*0.1)]);
    end

    % calculate distances
    parfor idx = 1:numel(linIdxs)
        % Convert back from a linear index to [i,j] co-ordinate pair
        [i,j] = ind2sub([nWv, nWv], linIdxs(idx));
        dist_mat(idx) = dist_fun(wv(i,:),wv(j,:));
    end
    
    % fix dist_mat shape from vector
    dist_mat = reshape(dist_mat,[nWv, nWv]);
end
% EOF