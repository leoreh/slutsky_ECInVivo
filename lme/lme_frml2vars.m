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
%   See also: LME_ANALYSE, REGEXP, UNIQUE
%
%   250106 Refactored for tidy separation of effects.
%   260108 Revised extraction logic to support complex formulas.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

% Ensure char format
frml = char(frml);


%% ========================================================================
%  EXTRACT RESPONSE VARIABLE
%  ========================================================================

% Get response variable (everything before the tilde ~)
% Matches start of string, valid var name, followed by ~
tokens = regexp(frml, '^\s*([a-zA-Z_]\w*)\s*~', 'tokens', 'once');

if isempty(tokens)
    error('lme_frml2vars:InvalidFormula', ...
        'Invalid formula format or missing response variable. Expected "Response ~ Predictors"');
end
varRsp = tokens{1};


%% ========================================================================
%  EXTRACT RANDOM EFFECTS
%  ========================================================================

% Regex to find terms grouped by parentheses containing a pipe |
% Matches: (ValidMatlabVar | ValidMatlabVar) or similar structures
regexRand = '\([^)]+\|[^)]+\)';
varsRand = regexp(frml, regexRand, 'match');

% Ensure column vector
if isempty(varsRand)
    varsRand = {};
else
    varsRand = varsRand(:)';
end


%% ========================================================================
%  EXTRACT FIXED EFFECTS
%  ========================================================================

% 1. Remove Response variable and Tilde from the formula
% (We use the regex derived from parsing to be precise)
rhs = regexprep(frml, '^\s*[a-zA-Z_]\w*\s*~', '');

% 2. Remove Random Effects terms
rhs = regexprep(rhs, regexRand, '');

% 3. Extract Variable Names
% Find all valid Matlab identifiers in the remaining string
% This handles A+B, A*B, A:B, etc.
rawVars = regexp(rhs, '[a-zA-Z_]\w*', 'match');

if isempty(rawVars)
    varsFxd = {};
else
    % Remove duplicates and sort
    varsFxd = unique(rawVars);

    % Ensure row vector
    varsFxd = varsFxd(:)';
end

end     % EOF