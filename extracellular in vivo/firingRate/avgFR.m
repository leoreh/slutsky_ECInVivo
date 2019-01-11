function avgfr = avgFR(spkcount, varargin)

% for each unit calculates the mean / max FR within a window
%
% INPUT
% required:
%   fr          struct (see calcFR)
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

if win(1) == 0; win(1) = 1; end
if win(2) == Inf; win(2) = size(spkcount.strd, 2); end

nunits = size(spkcount.strd, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

avgfr = zeros(nunits, 1);
switch method
    case 'max'
        avgfr = max(spkcount.strd(:, win(1) : win(2)), [], 2);
    case 'avg'
        avgfr = mean(spkcount.strd(:, win(1) : win(2)), 2);
end

end

% EOF