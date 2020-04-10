function x = remNear(varargin)

% given indices, only the first within dist samples will be kept.
% 
% INPUT
%   x           vector of indices 
%   dist        distance threshold
%   flip        logical. go last to first [1] or first to last [0]
% 
% OUTPUT
%   x           without near peaks
% 
% 22 nov 19 LH.
% 
% TO DO LIST
%   add option to determine dist occurding to percentage
%   apply for matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'x', []);
addOptional(p, 'dist', []);
addOptional(p, 'flip', false, @islogical);

parse(p, varargin{:})
x = p.Results.x;
dist = p.Results.dist;
flip = p.Results.flip;

x = x(:);   % make sure column vec

if flip
    x = flipud(x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j = 1;
while 1
        x([false(j, 1); abs(x(j + 1 : end) - x(j)) < dist]) = [];
        j = j + 1;        
    if j > length(x)
        return
    end
end

% EOF