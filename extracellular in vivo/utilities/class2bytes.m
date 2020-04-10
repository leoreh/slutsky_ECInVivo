function nbytes = class2bytes(x, varargin)

% returns the size in bytes of a data point of certain class.
% 
% INPUT:
%   x        any variable or specific string representing class.
%   var      logical. if true, x is treated as a variable and nbytes is
%            extracted from its size. if var is false, x is treated as a
%            class and nbytes is determined accordingly.
% 
% 10 apr 20 LH

p = inputParser;
addOptional(p, 'var', false, @islogical);

parse(p, varargin{:})
var = p.Results.var;

if var
    w = whos('x');
    nbytes = w.bytes / numel(x);
else
    switch x
        case {'uchar', 'unsigned char', 'schar', 'signed char', 'int8', 'integer*1', 'uint8'}
            nbytes = 1;
        case {'int16', 'integer*2', 'uint16'}
            nbytes = 2;
        case {'int32', 'integer*4', 'uint32', 'single', 'real*4', 'float32', 'real*4'}
            nbytes = 4;
        case {'int64', 'integer*8', 'uint64', 'double', 'real*8', 'float64', 'real*8'}
            nbytes = 8;
    end
end
end

% EOF