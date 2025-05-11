function frmlSimple = lme_frml2char(frml, varargin)

% simplifies linear mixed effects model formula by removing random effects
% and intercept terms
% e.g., 'Burst ~ 1 + Group + UnitType + (1|Mouse)' -> 'Burst~Group+UnitType'
%
% 10 jan 25 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'rm_rnd', true, @islogical);
addOptional(p, 'sfx', '', @ischar);
addOptional(p, 'pfx', '', @ischar);
addOptional(p, 'resNew', '', @ischar);

parse(p, varargin{:});
rm_rnd = p.Results.rm_rnd;
sfx = p.Results.sfx;
pfx = p.Results.pfx;
resNew = p.Results.resNew;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if formula is empty
if isempty(frml)
    error('Formula cannot be empty');
end

% Replace response variable if resNew is provided
if ~isempty(resNew)

    % Find the original response variable (text before '~')
    resOrig = regexp(frml, '^[^~]+', 'match', 'once');
    % Ensure the new response variable ends with a space if it doesn't already
    % and the part to replace also includes the space for cleaner replacement
    % or directly replace up to the tilde.
    frml = regexprep(frml, ['^' strtrim(resOrig) '\s*~'], [strtrim(resNew) ' ~']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simplify formula
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove random effects term (anything within parentheses)
if rm_rnd
    frmlSimple = regexprep(frml, '\s*\+ ([^)]*\)', '');
else
    frmlSimple = frml;
end

% remove intercept term (1 +) or (+ 1)
frmlSimple = regexprep(frmlSimple, '\s*1\s*\+\s*', '');
frmlSimple = regexprep(frmlSimple, '\s*\+\s*1\s*', '');

% remove all whitespace
% frmlSimple = regexprep(frmlSimple, '\s+', '');
frmlSimple = regexprep(frmlSimple, '\~', '~ ');

% remove trailing operators if they exist
frmlSimple = regexprep(frmlSimple, '[+\-*\/]\s*$', '');

if rm_rnd
    % replace * with x
    frmlSimple = regexprep(frmlSimple, '\*', ' x ');
else
    frmlSimple = regexprep(frmlSimple, '\*', ' * ');
end

% add suffix and prefix
frmlSimple = sprintf('%s%s%s', pfx, frmlSimple, sfx);

end

% EOF