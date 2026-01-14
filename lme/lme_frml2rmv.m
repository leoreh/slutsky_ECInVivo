function newFrml = lme_frml2rmv(frml, varRmv)
% LME_FRML2RMV Reconstructs a formula after removing a specific variable or term.
%
%   NEWFRML = LME_FRML2RMV(FRML, varRmv)
%   Removes a specific term from the LME formula.
%
%   LOGIC:
%       1. If varRmv is a Main Effect (e.g., 'A'):
%          - Remove 'A' everywhere.
%          - Remove any interaction checking 'A' (e.g., 'A:B').
%          - effectively: remove if contains 'A'.
%
%       2. If varRmv is an Interaction (e.g., 'A:B'):
%          - Remove ONLY the interaction term 'A:B'.
%          - Main effects 'A' and 'B' are preserved.
%
%   INPUTS:
%       frml      - (char) Original LME formula.
%       varRmv    - (char) Name of the variable or interaction to remove.
%                   Can be 'A', 'Group', 'A:B', or 'A:C'.
%
%   OUTPUTS:
%       newFrml   - (char) Reconstructed formula string.
%
%   See also: LME_ABLATION, LME_FRML2VARS

%% ========================================================================
%  PREPARATION
%  ========================================================================

frml = char(frml);
varRmv = char(varRmv);

% 1. Parse Formula Components
% We need to expand everything to explicit terms (A + B + A:B) to handle
% removal safely.
[~, varRsp, varsRand] = lme_frml2vars(frml);

% Get RHS (Fixed Effects)
rhs = regexprep(frml, '^\s*[a-zA-Z_]\w*\s*~', '');
rhs = regexprep(rhs, '\([^)]+\|[^)]+\)', ''); % Remove random effects

% Split into raw terms
% Expand terms (handling parentheses and interactions)
rawTerms = expand_frml(rhs);

%% ========================================================================
%  EXPAND TERMS (Convert A*B to A + B + A:B)
%  ========================================================================
expandedTerms = {};

for i = 1:numel(rawTerms)
    t = rawTerms{i};

    if contains(t, '*')
        % Expand A*B -> A, B, A:B
        atoms = strtrim(strsplit(t, '*'));

        % Add single atoms
        expandedTerms = [expandedTerms, atoms];

        % Add interaction (sorted A:B)
        if numel(atoms) > 1
            expandedTerms{end+1} = strjoin(sort(atoms), ':');
        end

    elseif contains(t, ':')
        % Already explicit interaction, just normalize sort
        atoms = strtrim(strsplit(t, ':'));
        expandedTerms{end+1} = strjoin(sort(atoms), ':');

    else
        % Main effect
        expandedTerms{end+1} = t;
    end
end

expandedTerms = unique(expandedTerms); % Remove duplicates

%% ========================================================================
%  FILTER TERMS
%  ========================================================================
finalTerms = {};

% Determine Removal Mode
isInteractionRmv = contains(varRmv, ':') || contains(varRmv, '*');

if isInteractionRmv
    % Mode: Remove Specific Interaction logic
    % Normalize varRmv to sorted colon format
    atomsRmv = strtrim(strsplit(varRmv, {':', '*'}));
    normRmv = strjoin(sort(atomsRmv), ':');

    for i = 1:numel(expandedTerms)
        t = expandedTerms{i};
        % Keep valid terms
        if ~strcmp(t, normRmv)
            finalTerms{end+1} = t;
        end
    end

else
    % Mode: Remove Main Effect logic
    % Remove the variable AND any interaction containing it
    for i = 1:numel(expandedTerms)
        t = expandedTerms{i};

        % Check if this term contains the variable
        % We split term into atoms to check exact match
        termAtoms = strtrim(strsplit(t, ':'));

        if ~ismember(varRmv, termAtoms)
            finalTerms{end+1} = t;
        end
    end
end

%% ========================================================================
%  RECONSTRUCT
%  ========================================================================

% Fixed Effects
if isempty(finalTerms)
    fixedPart = '1';
else
    fixedPart = strjoin(finalTerms, ' + ');
end

% Random Effects
if isempty(varsRand)
    randPart = '';
else
    randPart = [' + ' strjoin(varsRand, ' + ')];
end

newFrml = sprintf('%s ~ %s%s', varRsp, fixedPart, randPart);

end

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
