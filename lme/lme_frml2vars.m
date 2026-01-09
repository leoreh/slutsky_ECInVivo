function [varsFxd, varRsp, varsRand, varsIntr] = lme_frml2vars(frml)
% LME_FRML2VARS Extracts fixed effects, random effects, interactions and response.
%
%   [VARSFXD, VARRSP, VARSRAND, VARSINTR] = LME_FRML2VARS(FRML) parses a
%   Linear Mixed-Effects (LME) model formula string. It separates the
%   Response Variable (Y), the unique Fixed Effects atoms (X), the Random
%   Effects Terms (Z), and explicitly listed Interaction Terms.
%
%   INPUTS:
%       frml        - (char/string) LME formula (e.g., 'Y ~ A * B + (1|S)').
%
%   OUTPUTS:
%       varsFxd     - (cell) Unique Fixed Effects variables (atoms) (e.g., {'A', 'B'}).
%       varRsp      - (char) Name of the Response variable (e.g., 'Y').
%       varsRand    - (cell) Random Effects terms (e.g., {'(1|S)'}).
%       varsIntr    - (cell) Explicit interaction terms (e.g., {'A:B'}).
%
%   See also: LME_ANALYSE, REGEXP, UNIQUE
%
%   250109 Refactored to extract interaction terms for ablation.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

% Ensure char format
frml = char(frml);

%% ========================================================================
%  EXTRACT RESPONSE VARIABLE
%  ========================================================================

% Get response variable (everything before the tilde ~)
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
regexRand = '\([^)]+\|[^)]+\)';
varsRand = regexp(frml, regexRand, 'match');

% Ensure column vector
if isempty(varsRand)
    varsRand = {};
else
    varsRand = varsRand(:)';
end

%% ========================================================================
%  EXTRACT FIXED EFFECTS & INTERACTIONS
%  ========================================================================

% 1. Remove Response variable and random effects
rhs = regexprep(frml, '^\s*[a-zA-Z_]\w*\s*~', '');
rhs = regexprep(rhs, regexRand, '');

% 2. Split into terms by '+' to analyze structure
terms = strtrim(strsplit(rhs, '+'));
terms(strcmp(terms, '')) = [];

% 3. Process terms to find Atoms and Interactions
allAtoms = {};
varsIntr = {};

for i = 1:numel(terms)
    t = terms{i};

    if contains(t, '*')
        % "Product" interaction: A*B -> A + B + A:B
        % Extract atoms
        subAtoms = strtrim(strsplit(t, '*'));
        allAtoms = [allAtoms, subAtoms];

        % Generate all pairwise/higher-order interactions (simplified: just A:B for now)
        % For A*B*C, strict expansion is complex.
        % We will simplify checking: if '*' exists, it implies full factorial.
        % However, specifically for ablation, we usually care about the highest order
        % term or specific listed interactions.
        % Let's synthesize the standard interaction form 'A:B' for output.
        if numel(subAtoms) > 1
            % Create A:B string
            varsIntr{end+1} = strjoin(sort(subAtoms), ':');
        end

    elseif contains(t, ':')
        % Explicit interaction: A:B
        subAtoms = strtrim(strsplit(t, ':'));
        allAtoms = [allAtoms, subAtoms];

        % Store exact interaction term (sorted for consistency)
        varsIntr{end+1} = strjoin(sort(subAtoms), ':');

    else
        % Main effect
        allAtoms{end+1} = t;
    end
end

% 4. Unique Atoms (varsFxd)
varsFxd = unique(allAtoms);
varsFxd(strcmp(varsFxd, '1')) = []; % Remove intercept
if isempty(varsFxd), varsFxd = {}; else, varsFxd = varsFxd(:)'; end

% 5. Unique Interactions
if isempty(varsIntr)
    varsIntr = {};
else
    varsIntr = unique(varsIntr);
    varsIntr = varsIntr(:)';
end

end     % EOF