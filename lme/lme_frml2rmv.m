function newFrml = lme_frml2rmv(frml, varRmv)
% LME_FRML2RMV Reconstructs a formula after removing a specific variable or term.
%
%   NEWFRML = LME_FRML2RMV(FRML, varRmv)
%   Removes a specific term from both Fixed and Random effects in an LME formula.
%
%   LOGIC:
%       1. Fixed Effects:
%          - If varRmv is a Main Effect (e.g., 'A'): Removes 'A' and any
%            interaction containing 'A' (e.g., 'A:B').
%          - If varRmv is an Interaction (e.g., 'A:B'): Removes ONLY that
%            specific interaction term.
%
%       2. Random Effects:
%          - Iterates through each random term (e.g., '(A|G)').
%          - Removes varRmv from the design matrix (LHS) of the random term.
%          - If the Left-Hand Side (LHS) becomes empty, it implies an
%            Intercept-only random effect '(1|G)'.
%          - If varRmv matches the Grouping Variable (RHS), the entire
%            random term is removed.
%
%   INPUTS:
%       frml      - (char) Original LME formula.
%       varRmv    - (char) Name of the variable or interaction to remove.
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
% We use LME_FRML2VARS to robustly identify Random Effects.
[~, varRsp, varsRand] = lme_frml2vars(frml);

% 2. Extract Fixed Effects Part (RHS)
% Remove Response (LHS)
rhs = regexprep(frml, '^\s*[a-zA-Z_]\w*\s*~', '');
% Remove Random Effects strings (e.g., '(A|B)') to isolate Fixed Effects.
% Note: We assume standard parenthesis nesting for random effects.
rhs = regexprep(rhs, '\([^)]+\|[^)]+\)', '');

% 3. Identify Removal Mode
% If varRmv contains ':' or '*', we are in "Interaction Removal" mode.
isInteractionRmv = contains(varRmv, ':') || contains(varRmv, '*');


%% ========================================================================
%  PROCESS FIXED EFFECTS
%  ========================================================================

% 1. Expand Fixed Effects into initial terms (e.g. "A*B")
termsFixed = expand_frml(rhs);

% 2. Explode short-hand interactions (A*B -> A, B, A:B)
% This ensures we can safely target "A:B" without losing "A" or "B".
termsFixed = explode_interactions(termsFixed);

% 3. Filter Fixed Terms
termsFixed = filter_terms(termsFixed, varRmv, isInteractionRmv);


%% ========================================================================
%  PROCESS RANDOM EFFECTS
%  ========================================================================
% We need to inspect each random term (e.g., '(Slope|Group)') and remove
% the variable from the slope if present.

newVarsRand = {};

for i = 1:numel(varsRand)
    rTerm = varsRand{i}; % e.g. '(pBspk|Name)'

    % Strip outer parentheses
    cleanTerm = regexprep(rTerm, '^\(|\)$', '');

    % Split into Design (LHS) and Group (RHS)
    parts = strsplit(cleanTerm, '|');
    if numel(parts) < 2
        % Malformed term? Keep it as is to be safe.
        newVarsRand{end+1} = rTerm; %#ok<*AGROW>
        continue;
    end

    lhs = strtrim(parts{1}); % Slope / Design
    rhs = strtrim(parts{2}); % Group

    % Check 1: Is the variable the Grouping Variable?
    % If we are removing 'Name' and term is '(X|Name)', remove entire term.
    if strcmp(rhs, varRmv)
        continue; % Drop this random effect entirely
    end

    % Check 2: Remove variable from Slope (LHS)
    termsLhs = expand_frml(lhs);
    termsLhs = explode_interactions(termsLhs);
    termsLhsClean = filter_terms(termsLhs, varRmv, isInteractionRmv);

    % Reconstruct LHS
    if isempty(termsLhsClean)
        % If we removed everything (e.g., removing A from (A|G)),
        % assume it degrades to a Random Intercept (1|G).
        lhsNew = '1';
    else
        lhsNew = strjoin(unique(termsLhsClean), ' + ');
    end

    % Rebuild Term
    newVarsRand{end+1} = sprintf('(%s | %s)', lhsNew, rhs);
end


%% ========================================================================
%  RECONSTRUCT FORMULA
%  ========================================================================

% Fixed Part
if isempty(termsFixed)
    fixedPart = '1';
else
    fixedPart = strjoin(unique(termsFixed), ' + ');
end

% Random Part
if isempty(newVarsRand)
    randPart = '';
else
    randPart = [' + ' strjoin(newVarsRand, ' + ')];
end

newFrml = sprintf('%s ~ %s%s', varRsp, fixedPart, randPart);

end

%% ========================================================================
%  HELPER: TERM FILTERING
%  ========================================================================
function finalTerms = filter_terms(inputTerms, varRmv, isInteractionRmv)
% FILTER_TERMS Removes terms matching varRmv from the list.

finalTerms = {};

if isInteractionRmv
    % MODE: Remove Specific Interaction
    % Example: Remove 'A:B'. Keep 'A', 'B'. Remove ONLY 'A:B'.

    % Normalize varRmv to sorted colon format (e.g., 'B:A' -> 'A:B')
    atomsRmv = strtrim(strsplit(varRmv, {':', '*'}));
    normRmv = strjoin(sort(atomsRmv), ':');

    for i = 1:numel(inputTerms)
        t = inputTerms{i};

        % Normalize term to check match
        % e.g., 'B:A' -> 'A:B'
        if contains(t, ':') || contains(t, '*')
            atomsT = strtrim(strsplit(t, {':', '*'}));
            normT = strjoin(sort(atomsT), ':');
        else
            normT = t;
        end

        % Keep if NOT exact match
        if ~strcmp(normT, normRmv)
            finalTerms{end+1} = t;
        end
    end

else
    % MODE: Remove Main Effect
    % Example: Remove 'A'. Remove 'A' AND 'A:B'.

    for i = 1:numel(inputTerms)
        t = inputTerms{i};

        % Split term into atoms to find if it contains varRmv
        termAtoms = strtrim(strsplit(t, {':', '*'}));

        if ~ismember(varRmv, termAtoms)
            finalTerms{end+1} = t;
        end
    end
end

end

%% ========================================================================
%  HELPER: EXPLODE INTERACTIONS
%  ========================================================================
function newTerms = explode_interactions(terms)
% EXPLODE_INTERACTIONS Expands shorthand 'A*B' into 'A', 'B', 'A:B'.
% Also ensures plain interactions 'A:B' are normalized.

newTerms = {};
for i = 1:numel(terms)
    p = terms{i};

    if contains(p, '*')
        % A*B -> A, B, A:B
        subAtoms = strtrim(strsplit(p, '*'));
        newTerms = [newTerms, subAtoms]; %#ok<AGROW>

        if numel(subAtoms) > 1
            newTerms{end+1} = strjoin(sort(subAtoms), ':'); %#ok<AGROW>
        end
    elseif contains(p, ':')
        % Already explicit: A:B -> A:B (normalized)
        atoms = strtrim(strsplit(p, ':'));
        newTerms{end+1} = strjoin(sort(atoms), ':'); %#ok<AGROW>
    else
        % Main Effect
        newTerms{end+1} = p; %#ok<AGROW>
    end
end
newTerms = unique(newTerms);
end

%% ========================================================================
%  HELPER: FORMULA EXPANSION
%  ========================================================================
function terms = expand_frml(str)
% Recursively expands formula string into additive terms
% e.g. '(A+B)*C' -> {'A*C', 'B*C'} which eventually become {'A','B','A:B'}

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
                    % Concat with '*' to preserve operation for explode_interactions
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
% Checks if (A...B) is truly enclosed
if ~startsWith(str, '(') || ~endsWith(str, ')'), tf = false; return; end
count = 0;
for i = 1:length(str)-1
    if str(i) == '(', count = count + 1; end
    if str(i) == ')', count = count - 1; end
    if count == 0, tf = false; return; end
end
tf = true;
end
