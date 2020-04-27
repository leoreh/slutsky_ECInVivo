function [sig, dc] = rmDC(sig, varargin)

% removes dc offset from signal by calculating mean in window (baseline).
% adjusts class of dc to that of sig. function only works for 1/2
% dimensions.
% 
% INPUT
%   sig         numeric array. if matrix, baseline is calculated separatly
%               for each column (dim = 1) or row {dim = 2}
%   winCalc     time window for calculating DC {[1 Inf]}. specified as
%               index to mat rows.
%   dim         numeric. dimension on which to calculate mean {2}.
% 
% OUTPUT
%   sig         input singal without DC component   
%   dc          mean signal in window
% 
% 09 mar 19 LH
% 10 apr 20     compatibability with different classes and dim
% 
% notes
% 10 apr 20     compared performance on array of 35 x 5e6. if class double,
%               rmDC took 1.2 s but the remainder mean after dc removal was
%               <1e-8 - 1e-11. if single, rmDC took 0.83 s and the
%               remainder was <1e-2 - 1e-4. if int16, rmDc took 0.79 s and
%               the remainder was 0.1-0.5.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validate_win = @(win) assert(isnumeric(win) && length(win) == 2,...
    'time window must be in the format [start end]');

p = inputParser;
addOptional(p, 'win', [1 Inf], validate_win);
addOptional(p, 'dim', 2, @isnumeric);

parse(p,varargin{:})
win = p.Results.win;
dim = p.Results.dim;

if win(2) == inf
    win(2) = size(sig, dim);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc mean
if dim == 2
    dc = mean(sig(:, win(1) : win(2)), 2);
elseif dim == 1
    dc = mean(sig(win(1) : win(2), :), 1);
else
    error('dimension incorrect')
end
% convert to same class as sig
cl = class(sig);
eval(['dc = ' cl '(dc);'])
% remove mean from sig
sig = sig - dc;

end

% EOF