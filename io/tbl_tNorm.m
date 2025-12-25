function tblOut = tbl_tNorm(tbl, varsInc, winNorm, varsGrp)
% TBL_TNORM Normalizes vector columns in a table to a specific window.
%
% SUMMARY:
% This function takes a table containing vector columns (either as matrices
% or cell arrays of vectors) and normalizes them. Normalization is performed
% by dividing the vector by the mean value within a specified window index.
% This can be done row-by-row (each row normalized to its own baseline) or
% group-wise (each row normalized to the group's grand mean baseline).
%
% INPUT:
%   tbl         - Input table.
%   varsInc     - Cell array of variable names (strings) to normalize.
%                  Columns must be numeric matrices (N x M) or cell arrays
%                  of numeric vectors.
%   winNorm     - 2-element numeric vector [startIndex, endIndex] defining
%                  the normalization window (inclusive, 1-based indices).
%   varsGrp     - (Optional) Cell array of categorical variable names for
%                  grouping.
%                  * If EMPTY (default): Row-wise normalization. Each row
%                    is divided by mean(rowVector(winNorm)).
%                  * If SPECIFIED: Group-wise normalization. Each row in a
%                    group is divided by the GRAND MEAN of all vectors in
%                    that group within the window.
%
% OUTPUT:
%   tblOut      - Table with normalized columns.
%
% EXAMPLE:
%   % Row-wise normalization of 'LFP' to first 100 points
%   tbl = tbl_tNorm(tbl, {'LFP'}, [1, 100], {});
%
%   % Group-wise normalization of 'FiringRate' by 'Genotype'
%   tbl = tbl_tNorm(tbl, {'FiringRate'}, [1, 300], {'Genotype'});
%
% DEPENDENCIES:
%   None
%
%   See also: TBL_TRANSFORM

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

if nargin < 4
    varsGrp = {};
end

% Input validation
if ~istable(tbl)
    error('Input "tbl" must be a table.');
end

if isstring(varsInc)
    varsInc = cellstr(varsInc);
end
if ~iscell(varsInc)
    error('"varsInc" must be a cell array of strings.');
end

if ~isnumeric(winNorm) || numel(winNorm) ~= 2
    error('"winNorm" must be a 2-element numeric vector [start, end].');
end

if isstring(varsGrp)
    varsGrp = cellstr(varsGrp);
end
if ~iscell(varsGrp)
    error('"varsGrp" must be a cell array of strings (or empty for row-wise).');
end

tblOut = tbl;

% Start/End indices
idxStart = winNorm(1);
idxEnd = winNorm(2);


%% ========================================================================
%  GROUP INDICES
%  ========================================================================

if ~isempty(varsGrp)
    % Find unique combinations of grouping variables
    try
        uGrps = unique(tblOut(:, varsGrp), 'rows');
    catch ME
        error('Failed to group by varsGrp. Ensure columns exist and are categorical/numeric. Msg: %s', ME.message);
    end

    % Map each row to a group index
    [~, loc] = ismember(tblOut(:, varsGrp), uGrps, 'rows');
    % loc contains index of uGrps for each row in tbl
else
    % No grouping -> Row-wise operations don't need group indices in the same way,
    % but we handle logic differently below.
    loc = [];
end


%% ========================================================================
%  TRANSFORM
%  ========================================================================

for iVar = 1:length(varsInc)
    varName = varsInc{iVar};

    if ~ismember(varName, tbl.Properties.VariableNames)
        warning('Variable "%s" not found in table. Skipping.', varName);
        continue;
    end

    dataRaw = tblOut.(varName);

    % Determine data type (Matrix or Cell Array)
    isCell = iscell(dataRaw);

    if ~isCell && ~isnumeric(dataRaw)
        warning('Variable "%s" is neither numeric matrix nor cell array. Skipping.', varName);
        continue;
    end

    % ---------------------------------------------------------
    % LOGIC BRANCH 1: ROW-WISE NORMALIZATION (varsGrp is empty)
    % ---------------------------------------------------------
    if isempty(varsGrp)

        if isCell
            % Cell Array of Vectors
            dataNew = dataRaw;
            for iRow = 1:height(tblOut)
                vec = dataRaw{iRow};
                if isempty(vec) || length(vec) < idxEnd
                    % Handle short/empty vectors if necessary
                    continue;
                end

                blVal = mean(vec(idxStart:idxEnd), 'all', 'omitnan');

                % Avoid division by zero
                if blVal == 0, blVal = eps; end

                dataNew{iRow} = vec ./ blVal;
            end
            tblOut.(varName) = dataNew;

        else
            % Numeric Matrix (N x M)
            % Check bounds
            if size(dataRaw, 2) < idxEnd
                warning('Variable "%s" width (%d) is smaller than window end (%d).', ...
                    varName, size(dataRaw, 2), idxEnd);
            end

            % Calculate baseline for each row
            % winNorm indices apply to columns
            colsIdx = max(1, idxStart) : min(size(dataRaw, 2), idxEnd);
            blVals = mean(dataRaw(:, colsIdx), 2, 'omitnan');

            % Handle zeros
            blVals(blVals == 0) = eps;

            % Divide
            tblOut.(varName) = dataRaw ./ blVals;
        end

        % ---------------------------------------------------------
        % LOGIC BRANCH 2: GROUP-WISE NORMALIZATION
        % ---------------------------------------------------------
    else

        dataNew = dataRaw; % Initialize output structure

        nGrps = height(uGrps);

        for iGrp = 1:nGrps
            % Rows belonging to this group
            idxRows = (loc == iGrp);

            if ~any(idxRows), continue; end

            % Extract all baseline data for this group to calculate ONE Grand Mean
            allBlValues = [];

            if isCell
                % Extract relevant window from each cell in the group
                groupCells = dataRaw(idxRows);
                for k = 1:length(groupCells)
                    vec = groupCells{k};
                    if ~isempty(vec) && length(vec) >= idxEnd
                        allBlValues = [allBlValues, vec(idxStart:idxEnd)]; %#ok<AGROW>
                    elseif ~isempty(vec) && length(vec) >= idxStart
                        % Partial overlap
                        allBlValues = [allBlValues, vec(idxStart:end)]; %#ok<AGROW>
                    end
                end
            else
                % Matrix: Extract submatrix
                groupMat = dataRaw(idxRows, :);
                colsIdx = max(1, idxStart) : min(size(groupMat, 2), idxEnd);
                if ~isempty(colsIdx)
                    % Flatten the windowed part of the matrix
                    allBlValues = groupMat(:, colsIdx);
                    allBlValues = allBlValues(:);
                end
            end

            % Calculate Grand Mean
            grandMean = mean(allBlValues, 'all', 'omitnan');

            if isempty(grandMean) || isnan(grandMean) || grandMean == 0
                grandMean = eps;
            end

            % Normalize rows in this group by the Grand Mean
            if isCell
                % Apply to each cell individually
                currCells = dataRaw(idxRows);
                for k = 1:length(currCells)
                    if ~isempty(currCells{k})
                        currCells{k} = currCells{k} ./ grandMean;
                    end
                end
                dataNew(idxRows) = currCells;
            else
                % Matrix division
                dataNew(idxRows, :) = dataRaw(idxRows, :) ./ grandMean;
            end

        end

        tblOut.(varName) = dataNew;

    end

end

end     % EOF
