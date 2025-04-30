function [lme_tbl, norm_cfg] = mcu_lmeNorm(lme_tbl, norm_var, norm_val)
% mcu_lmeNorm normalizes measurement data in an LME table relative to a baseline.
%
% This function takes an lme_tbl (typically produced by mcu_lmeOrg) and
% normalizes its measurement variable (the first variable in the table)
% for each mouse by dividing each data point by the mean value when the
% variable specified by norm_var equals norm_val. For example, if norm_var
% is 'Day' and norm_val is 'BSL', then each measurement for a given mouse
% is divided by the mouse's mean measurement at baseline.
%
% INPUT
%   lme_tbl   Table organized for LME analysis. The first variable is assumed
%             to be the measurement (e.g., FOOOF parameter) and the table must
%             contain a column with the name given in norm_var (e.g., 'Day').
%   norm_var  Char specifying the variable in lme_tbl to use for baseline grouping.
%             For example, 'Day'.
%   norm_val  Char or categorical value specifying the baseline value in norm_var.
%             For example, 'BSL'.
%
% OUTPUT
%   lme_tbl   Table with normalized measurement values.
%   norm_cfg  Struct containing normalization configuration including:
%             norm_var, norm_val, and a table of baseline means per mouse.
%
% 07 Feb 25

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    error('mcu_lmeNorm:NotEnoughInputs', 'Function requires lme_tbl, norm_var, and norm_val.');
end

if ~ismember(norm_var, lme_tbl.Properties.VariableNames)
    error('mcu_lmeNorm:MissingVar', 'The variable "%s" is not present in the table.', norm_var);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identify measurement variable and unique mice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assume the first variable in the table holds the measurement values.
yVar = lme_tbl.Properties.VariableNames{1};

% Ensure the grouping variable 'Mouse' exists.
if ~ismember('Mouse', lme_tbl.Properties.VariableNames)
    error('mcu_lmeNorm:MissingMouse', 'The table must contain a "Mouse" variable.');
end

% Get unique mouse identifiers.
mice = unique(lme_tbl.Mouse);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalize measurement data per mouse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate vector to store each mouse's baseline mean.
baselineMeanVec = nan(size(mice));

for i = 1:length(mice)
    curr_mouse = mice(i);
    
    % Find indices for current mouse.
    idx_mouse = lme_tbl.Mouse == curr_mouse;
    
    % Find indices for baseline condition (e.g., Day equals 'BSL').
    idx_baseline = idx_mouse & (lme_tbl.(norm_var) == norm_val);
    
    % Compute baseline mean for the measurement variable.
    if any(idx_baseline)
        baseMean = mean(lme_tbl.(yVar)(idx_baseline) + eps);
    else
        warning('mcu_lmeNorm:NoBaseline', ...
            'Mouse %s has no entries where %s == %s. Normalization skipped for this mouse.', ...
            string(curr_mouse), norm_var, norm_val);
        baseMean = NaN;
    end
    baselineMeanVec(i) = baseMean;
    
    % Normalize all measurements for this mouse if a valid baseline is available.
    if ~isnan(baseMean) && baseMean ~= 0
        lme_tbl.(yVar)(idx_mouse) = lme_tbl.(yVar)(idx_mouse) / baseMean;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store normalization configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

norm_cfg.norm_var = norm_var;
norm_cfg.norm_val = norm_val;
norm_cfg.baselineMeans = table(mice, baselineMeanVec, 'VariableNames', {'Mouse', 'BaselineMean'});

end
% EOF
