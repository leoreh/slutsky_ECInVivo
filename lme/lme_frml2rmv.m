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
rawTerms = strtrim(strsplit(rhs, '+'));
rawTerms(strcmp(rawTerms, '')) = [];

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
