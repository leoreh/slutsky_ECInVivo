function avgfr = avgFR(fr, varargin)

% for each unit calculates the mean / max FR within a window
%
% INPUT
% required:
%   fr          matrix of units (rows) x firing rate in time bins (columns)
%               for example, fr.strd from calcFR is a valid input.
% optional:
%   win         time window for calculation {[1 Inf]}. specified in min.
%   method      calculate 'max' or 'avg' FR within win {'avg'}.
%
% OUTPUT
% avgfr         vector with mean / max FR across units (rows)
%
% 11 jan 19 LH.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validate_win = @(win) assert(isnumeric(win) && length(win) == 2,...
    'time window must be in the format [start end]');

p = inputParser;
addOptional(p, 'win', [1 Inf], validate_win);
addOptional(p, 'method', 'avg', @ischar);

parse(p,varargin{:})
win = p.Results.win;
method = p.Results.method;

[nunits, nmints] = size(fr);
if win(1) == 0; win(1) = 1; end
if win(2) == Inf; win(2) = nmints; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
avgfr = zeros(nunits, 1);
switch method
    case 'max'
        avgfr = max(fr(:, win(1) : win(2)), [], 2);
    case 'avg'
        avgfr = mean(fr(:, win(1) : win(2)), 2);
end

end

% EOF