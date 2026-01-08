function newFrml = lme_frml2rmv(frml, varRmv)
% LME_FRML2RMV Reconstructs a formula after removing a specific variable.
%
%   NEWFRML = LME_FRML2RMV(FRML, varRmv)
%    parses the input formula FRML, removes the specified variable varRmv
%   from all terms (handling interactions and nested structures), and returns
%   the reconstructed formula string NEWFRML.
%
%   This is critical for ablation analysis where removing a main effect (e.g.,
%   'Group') from an interaction term ('fr * Group') should result in the
%   retention of the partner variable ('fr'), rather than removing the entire
%   interaction term or leaving a dangling operator.
%
%   LOGIC:
%       1. Interaction Terms (A * B):
%          - Removing A leaves B.
%          - Removing B leaves A.
%       2. Explicit Interactions (A : B):
%          - Removing A removes the entire term (A:B requires A).
%       3. Main Effects (A):
%          - Removing A removes the term.
%
%   INPUTS:
%       frml        - (char) Original LME formula.
%       varRmv    - (char) Name of the variable to remove.
%
%   OUTPUTS:
%       newFrml     - (char) New formula with the variable removed.
%
%   EXAMPLES:
%       f = 'y ~ A * B + (1|S)';
%       lme_frml2rmv(f, 'A') -> 'y ~ B + (1|S)'
%       lme_frml2rmv(f, 'B') -> 'y ~ A + (1|S)'
%
%       f = 'y ~ A + A:B + (1|S)';
%       lme_frml2rmv(f, 'A') -> 'y ~ (1|S)'  (A:B depends on A)
%
%   See also: LME_ABLATION, LME_FRML2VARS

%% ========================================================================
%  PREPARATION
%  ========================================================================

frml = char(frml);
varRmv = char(varRmv);

% 1. Separate Random Effects and Response
[~, varRsp, varsRand] = lme_frml2vars(frml);

% 2. Extract the Fixed Effects string (between ~ and random effects)
% Use lme_frml2vars logic to isolate the RHS, then strip random effects
rhs = regexprep(frml, '^\s*[a-zA-Z_]\w*\s*~', '');
regexRand = '\([^)]+\|[^)]+\)';
fixedStr = regexprep(rhs, regexRand, '');

% 3. Split by '+' to get individual terms
% Handle potential whitespace
terms = strtrim(strsplit(fixedStr, '+'));
terms(strcmp(terms, '')) = []; % Remove empty


%% ========================================================================
%  PROCESS TERMS
%  ========================================================================

newTerms = {};

for iVar = 1:length(terms)
    term = terms{iVar};

    % Check for Interaction (*)
    if contains(term, '*')
        % "Product" interaction: A * B -> A + B + A:B
        % Rule: Remove varRmv from the list of atoms.
        % If A * B * C and we remove B, we get A * C.
        atoms = strtrim(strsplit(term, '*'));
        atoms(strcmp(atoms, varRmv)) = [];

        if ~isempty(atoms)
            newTerms{end+1} = strjoin(atoms, ' * ');
        end

        % Check for Explicit Interaction (:)
    elseif contains(term, ':')
        % "Colon" interaction: A : B
        % Rule: If varRmv is present, the specific interaction is invalid.
        atoms = strtrim(strsplit(term, ':'));
        if ~any(strcmp(atoms, varRmv))
            newTerms{end+1} = term;
        end

        % Main Effect / Single Variable
    else
        if ~strcmp(term, varRmv)
            newTerms{end+1} = term;
        end
    end
end


%% ========================================================================
%  RECONSTRUCT FORMULA
%  ========================================================================

% 1. Fixed Effects Part
if isempty(newTerms)
    fixedPart = '1';
else
    fixedPart = strjoin(newTerms, ' + ');
end

% 2. Random Effects Part
if isempty(varsRand)
    randPart = '';
else
    randPart = [' + ' strjoin(varsRand, ' + ')];
end

% 3. Combine
newFrml = sprintf('%s ~ %s%s', varRsp, fixedPart, randPart);

end
