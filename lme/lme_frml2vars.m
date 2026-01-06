function [varsFxd, varRsp, varsRand] = lme_frml2vars(frml)
% LME_FRML2VARS Extracts fixed effects, random effects, and response variable.
%
%   [VARSFXD, VARRSP, VARSRAND] = LME_FRML2VARS(FRML) parses a Linear
%   Mixed-Effects (LME) model formula string. It separates the Response
%   Variable (Y), the Fixed Effects Predictors (X), and the Random Effects
%   Terms (Z).
%
%   INPUTS:
%       frml        - (char/string) LME formula (e.g., 'Y ~ A * B + (1|Subject)').
%
%   OUTPUTS:
%       varsFxd     - (cell) Names of Fixed Effects variables (e.g., {'A', 'B'}).
%       varRsp      - (char) Name of the Response variable (e.g., 'Y').
%       varsRand    - (cell) Random Effects terms (e.g., {'(1|Subject)'}).
%
%   EXAMPLES:
%       [vf, vr, vz] = lme_frml2vars('Spikes ~ Depth + (1|Mouse)');
%       % vf = {'Depth'}, vr = 'Spikes', vz = {'(1|Mouse)'}
%
%   See also: REGEXP, STRTRIM, UNIQUE
%
%   250106 Refactored for tidy separation of effects.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

% Ensure char format
frml = char(frml);

%% ========================================================================
%  EXTRACT RESPONSE VARIABLE
%  ========================================================================

% Get response variable (everything before the tilde ~)
varRsp = regexp(frml, '^\s*(\w+)\s*~', 'tokens', 'once');

if isempty(varRsp)
    error('lme_frml2vars:InvalidFormula', ...
        'Invalid formula format. Expected "Response ~ Predictors"');
end
varRsp = varRsp{1};


%% ========================================================================
%  EXTRACT RANDOM EFFECTS
%  ========================================================================

% Regex to find terms grouped by parentheses containing a pipe |
% Matches: (ValidMatlabVar | ValidMatlabVar) or similar structures
regexRand = '\([^)]+\|[^)]+\)';
varsRand = regexp(frml, regexRand, 'match');

% Remove random effects from the formula string to isolate fixed effects
xStr = regexprep(frml, ['^\s*\w+\s*~|' regexRand], '');


%% ========================================================================
%  EXTRACT FIXED EFFECTS
%  ========================================================================

% Initialize
varsFxd = {};

if isempty(strtrim(xStr))
    % Case: Only Intercept or Empty Fixed Effects
    return;
end

% Check for interaction terms (A * B format)
interactVars = regexp(xStr, '(\w+)\s*[*]\s*(\w+)', 'tokens');

if ~isempty(interactVars)
    % Process all interaction pairs
    for iVar = 1:length(interactVars)
        varsFxd = [varsFxd, interactVars{iVar}];
    end
else
    % No interactions: standard additive model

    % Split by '+'
    xVars = strtrim(strsplit(xStr, '+'));

    % Remove '1' (Intercept)
    xVars(strcmp(xVars, '1')) = [];

    % Filter empty cells
    varsFxd = xVars(~cellfun(@isempty, xVars));
end

% Remove duplicates (e.g. from A*B -> A and B) and sort
varsFxd = unique(varsFxd, 'stable');

end     % EOF