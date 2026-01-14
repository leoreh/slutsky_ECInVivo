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
% 2. Expand terms (handling parentheses and interactions)
terms = expand_frml(rhs);

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

        % Generate standard interaction form 'A:B'
        if numel(subAtoms) > 1
            varsIntr{end+1} = strjoin(sort(subAtoms), ':');
        end

    elseif contains(t, ':')
        % Explicit interaction: A:B
        % Do NOT add to allAtoms (Fixed Effects).
        % Because A:B does not imply A or B is a Main Effect.
        subAtoms = strtrim(strsplit(t, ':'));

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

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function terms = expand_frml(str)
% Recursively expands formula string into additive terms
% e.g. '(A+B)*C' -> {'A*C', 'B*C'}

str = strtrim(str);
if isempty(str), terms = {}; return; end

% 1. Split by top-level '+'
parts = split_balanced(str, '+');

terms = {};
for i = 1:numel(parts)
    p = strtrim(parts{i});
    if isempty(p), continue; end

    % 2. Handle Multiplication/Distribution
    factors = split_balanced(p, '*');

    if numel(factors) == 1
        % No multiplication at this level
        % Check for surrounding parens: (A+B)
        if startsWith(p, '(') && endsWith(p, ')') && is_enclosed(p)
            % Strip parens and recurse
            inner = p(2:end-1);
            terms = [terms, expand_frml(inner)]; %#ok<AGROW>
        else
            % Atom or Interaction (A:B) or plain var
            terms{end+1} = p; %#ok<AGROW>
        end
    else
        % Distribute factors: (A+B)*C -> expand(A+B) X expand(C)
        expandedFactors = cell(size(factors));
        for k = 1:numel(factors)
            expandedFactors{k} = expand_frml(factors{k});
        end

        % Cartesian product of expanded factors
        prodTerms = expandedFactors{1};
        for k = 2:numel(factors)
            nextSet = expandedFactors{k};
            newProd = {};
            for x = 1:numel(prodTerms)
                for y = 1:numel(nextSet)
                    newProd{end+1} = [prodTerms{x} '*' nextSet{y}]; %#ok<AGROW>
                end
            end
            prodTerms = newProd;
        end
        terms = [terms, prodTerms]; %#ok<AGROW>
    end
end
end

function parts = split_balanced(str, delim)
% Splits string by delimiter, respecting parentheses
parts = {};
current = '';
parenLvl = 0;

for i = 1:length(str)
    c = str(i);
    if c == '('
        parenLvl = parenLvl + 1;
    elseif c == ')'
        parenLvl = parenLvl - 1;
    end

    if c == delim && parenLvl == 0
        parts{end+1} = current; %#ok<AGROW>
        current = '';
    else
        current = [current c]; %#ok<AGROW>
    end
end
parts{end+1} = current;
end

function tf = is_enclosed(str)
% Checks if (A...B) is truly enclosed or just A starts with ( and B ends with )
% e.g. (A)+(B) -> False
%      (A+B)   -> True
if ~startsWith(str, '(') || ~endsWith(str, ')'), tf = false; return; end
count = 0;
for i = 1:length(str)-1
    if str(i) == '(', count = count + 1; end
    if str(i) == ')', count = count - 1; end
    if count == 0, tf = false; return; end % Closed before end
end
tf = true;
end