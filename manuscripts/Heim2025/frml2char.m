function frmlSimple = frml2char(frml, varargin)

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

parse(p, varargin{:});
rm_rnd = p.Results.rm_rnd;
sfx = p.Results.sfx;
pfx = p.Results.pfx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if formula is empty
if isempty(frml)
    error('Formula cannot be empty');
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