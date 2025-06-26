function [varsFxd, varRsp] = lme_frml2vars(frml, varargin)
% LME_FRML2VARS Extracts fixed effects and response variable names from LME formula.
%
% This function parses a linear mixed-effects model formula string to extract
% the response variable and fixed effects variables, handling interaction terms
% and removing random effects and intercept terms.
%
% INPUTS:
%   frml        (Required) Char/String: LME formula (e.g., 'Y ~ A * B + (1|Subject)').
%
% VARARGIN (Name-Value Pairs):
%   'flgRand'   (Optional) Logical: Whether to include intercept term.
%                Default is false.
%
% OUTPUTS:
%   varsFxd     Cell array: Names of fixed effects variables.
%   varRsp      Char: Name of the response variable.
%
% EXAMPLES:
%   [varsFxd, varRsp] = lme_frml2vars('Y ~ A * B + (1|Subject)')
%   % Returns: varsFxd = {'A', 'B'}, varRsp = 'Y'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p, 'flgRand', false, @islogical);
parse(p, varargin{:});

opts = p.Results;
frml = char(frml); % Ensure char format

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACT RESPONSE VARIABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get response variable (everything before ~)
varRsp = regexp(frml, '(\w+)\s*~', 'tokens');
if isempty(varRsp)
    error('Invalid formula format. Expected format: "Response ~ Predictors"');
end
varRsp = varRsp{1}{1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACT FIXED EFFECTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get everything between ~ and ( or end of string
xStr = regexp(frml, '~\s*(.*?)\s*(\(|$)', 'tokens');
if isempty(xStr)
    error('Invalid formula format. No fixed effects found.');
end
xStr = xStr{1}{1}; % Extract matched string

% Initialize varsFxd as empty cell array
varsFxd = {};

% First check for interaction terms (A * B format)
interactVars = regexp(xStr, '(\w+)\s*[*]\s*(\w+)', 'tokens');
if ~isempty(interactVars)
    % Process all interaction pairs
    for iVar = 1:length(interactVars)
        varsFxd = [varsFxd, interactVars{iVar}];
    end
    % Remove duplicates while preserving order
    varsFxd = unique(varsFxd, 'stable');
else
    % If no interactions, process normally
    if ~opts.flgRand
        xStr = regexprep(xStr, '\s*1\s*\+?\s*', ''); % Remove intercept term
    end
    xStr = regexprep(xStr, '\([^)]*\)', ''); % Remove random effects
    xVars = strtrim(strsplit(xStr, '+')); % Split by + and trim whitespace
    varsFxd = xVars(~cellfun(@isempty, xVars)); % Remove any empty cells
end

end 