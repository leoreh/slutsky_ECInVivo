function [frnorm, ithr, istable] = normFR(fr, varargin)

% normalizes FR according to maximum / average in a certain window.
% applies thershold to select high firing and/or stable units. 
% 
% INPUT
% required:
%   fr          matrix of units (rows) x firing rate in time bins (columns)
% optional:
%   metBL       method to calculate baseline as 'max' or {'avg'}.
%   win         time window for calculation {[]}. index to fr.
%   select      cell array with strings expressing method to select units.
%               'thr' - units with fr > 0.05 during baseline
%               'stable' - units with fr std < fr avg during baseline. 
%               default = none.
%
% OUTPUT        
%   frnorm      matrix of units (rows) x normalized firing rate in
%               time bins (columns).
%   ithr        logical vector of units with FR above thr
%   istable     logical vector of units with stable firing
%
% 26 feb 19 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validate_win = @(win) assert(isnumeric(win) && length(win) == 2,...
    'time window must be in the format [start end]');

p = inputParser;
addParameter(p, 'metBL', 'avg', @ischar);
addParameter(p, 'win', [], validate_win);
addParameter(p, 'select', {}, @iscell);

parse(p, varargin{:})
metBL = p.Results.metBL;
win = p.Results.win;
select = p.Results.select;

[nunits, nbins] = size(fr);
if isempty(win); win = [1 nbins]; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalize FR to baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate baseline
bl = avgFR(fr, 'method', metBL, 'win', win);

% select units who fired above thr
if any(strcmp(select, 'thr'))
    ithr = bl > 0.05;
else
    ithr = ones(nunits, 1);
end

% select units with low variability
if any(strcmp(select, 'stable'))
    bl_std = std(fr(:, win(1) : win(2)), [], 2);
    istable = bl_std < bl;
else
    istable = ones(nunits, 1);
end

idx = istable & ithr;

% normalize
frnorm = fr(idx, :) ./ bl(idx);

end

% EOF