function norm_tbl = lme_normalize(varargin)
% LME_NORMALIZE Normalizes LME table data relative to reference categories.
%
% SUMMARY:
% This function normalizes the response variable in an LME table (lme_tbl)
% relative to the mean value of the reference category specified by normVar.
% The normalization is performed separately for each unique combination of
% categories defined in groupVars.
%
% INPUT:
%   lme_tbl     (Required) Table output from lme_org.
%   normVar     (Required) String. The name of the categorical variable whose
%               reference category mean will be used for normalization (e.g., 'Day').
%   groupVars   (Required) Cell array of strings. Names of categorical
%               variables defining the groups within which normalization
%               should occur (e.g., {'Group', 'State'}).
%
% OUTPUT:
%   norm_tbl    Table. The input table with the response variable normalized.
%               The normalized column replaces the original response variable column.
%
% EXAMPLE:
%   % Normalize FR by the mean FR during 'BSL' day, separately for each 'Group' and 'State'
%   normTable = lme_normalize('lme_tbl', myLmeTable, 'normVar', 'Day', 'groupVars', {'Group', 'State'});
%
% HISTORY:
%   07 Aug 24 - Initial version

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define input parameters and parse arguments
p = inputParser;
addOptional(p, 'lme_tbl', table()); % Use table() as a placeholder default
addOptional(p, 'normVar', '');
addOptional(p, 'groupVars', {});
addOptional(p, 'lme_cfg', struct()); % Add lme_cfg
parse(p, varargin{:});

% Extract parameters
lme_tbl = p.Results.lme_tbl;
normVar = p.Results.normVar;
groupVars = p.Results.groupVars;
lme_cfg = p.Results.lme_cfg; % Extract lme_cfg (though not used in core logic below)

% Ensure groupVars is cellstr
if isstring(groupVars)
    groupVars = cellstr(groupVars);
end

% Identify response variable (assumed to be the first column)
yName = lme_tbl.Properties.VariableNames{1};

% Get reference category (assuming first category is reference)
refCategory = categories(lme_tbl.(normVar));
refCategory = refCategory{1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NORMALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

norm_tbl = lme_tbl; % Copy table to modify
norm_data = nan(height(norm_tbl), 1); % Initialize normalized data vector

% Find unique combinations of grouping variables
uniqueGroups = unique(norm_tbl(:, groupVars), 'rows');

% Loop through each unique group combination
for i = 1:height(uniqueGroups)
    uRow = uniqueGroups(i, :);

    % Find rows matching the current unique group combination
    idx_group = true(height(norm_tbl), 1);
    for j = 1:length(groupVars)
        idx_group = idx_group & (norm_tbl.(groupVars{j}) == uRow.(groupVars{j}));
    end

    % Find rows within this group that match the reference category of normVar
    idx_ref = norm_tbl.(normVar) == refCategory;
    idx_ref_group = idx_group & idx_ref;

    % Calculate the mean for the reference category within this group
    refMean = mean(norm_tbl.(yName)(idx_ref_group), 'omitnan');

    % Check for issues with reference mean
    if isnan(refMean) || refMean == 0
        if refMean == 0
           refMean = eps; % Avoid division by zero, result will be large
        end
    end

    % Normalize the data for the current group (all normVar categories
    % within the group)
    norm_data(idx_group) = norm_tbl.(yName)(idx_group) / abs(refMean);

    % Set the reference group to exactly 1 (so after scaling it's 100%)
    norm_data(idx_ref_group) = 1;
end

% Replace the original response variable column with the normalized data
norm_tbl.(yName) = norm_data * 100;

% Optional: Rename the column to indicate normalization
% norm_tbl.Properties.VariableNames{yName} = [yName ' (Norm.)'];

end % end function lme_normalize
