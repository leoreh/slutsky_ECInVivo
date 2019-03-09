function mat = times2binary(timestamps, varargin)

% OBSOLETE! replaced by calcFR (w/ small binsize)
% converts a vector of timestamps to a continuous binary vector.
%
% INPUT
%   timestamps  a cell array of vectors. each vector (unit) contains the
%               timestamps of spikes. for example {spikes.times{1:4}}
%   win         only spikes within the window will be used.
%
% 14 jan 19 LH.
% 
% TO DO LIST
%   answer for cases where multiple spikes fall in the same bin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validate_win = @(win) assert(isnumeric(win) && length(win) == 2,...
    'time window must be in the format [start end]');

p = inputParser;
addOptional(p, 'win', [21600 21650], validate_win);

parse(p, varargin{:})
win = p.Results.win;

% adjust window
if win(2) == Inf
    for i = 1 : length(timestamps)
        recDur(i) = max(timestamps{i}(:, 1));
    end
    win(2) = max(recDur);
end

% constants
rfactor = 1e3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
nunits = length(timestamps);
mat = zeros(nunits, diff(win) * rfactor);
win = win * rfactor;

% convert
for i = 1 : nunits
    stamps = ceil(timestamps{i} * rfactor);
    if length(unique(stamps)) ~= length(timestamps{i})
        warning(' clu%d: %d spikes fall in the same timebin', i, length(timestamps{i}) - length(unique(stamps)))
    end
    idx = stamps(stamps > win(1) & stamps < win(end)) - win(1);
    for j = 1 : length(idx)
        x(j) = sum(stamps == idx(j))
    end
    mat(i, idx) = 1;   
    % mat(i, idx) = i; % useful for raster plot in prism 
end

% add time vector
mat = [0 : 1 / rfactor : diff(win) / rfactor - 1 / rfactor; mat]';

end

% EOF

