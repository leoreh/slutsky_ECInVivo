function [sig, dc] = rmDC(sig, varargin)

% removes DC from signal by calculating mean in window (baseline)
%  
% INPUT
%   sig         numeric array. if matrix, baseline is calculated separatly
%               for each column
%   winCalc     time window for calculating DC {[1 Inf]}. specified as
%               index to mat rows.
% 
% OUTPUT
%   sig         input singal without DC component   
%   dc          mean signal in window
% 
% 09 mar 19 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validate_win = @(win) assert(isnumeric(win) && length(win) == 2,...
    'time window must be in the format [start end]');

p = inputParser;
addOptional(p, 'win', [1 Inf], validate_win);

parse(p,varargin{:})
win = p.Results.win;

if win(2) == inf
    win(2) = size(sig, 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dc = mean(sig(win(1) : win(2), :));
sig = sig - dc;

end