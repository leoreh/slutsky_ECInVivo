function varargout = symlog(varargin)
% SYMLOG Symmetric Logarithmic Transformation (Bi-Symmetric Log)
%
%   Y = SYMLOG(X) transforms the numeric array X using the symmetric
%   logarithmic transformation. This is useful for data that spans several
%   orders of magnitude and includes both positive and negative values (and
%   zero).
%
%       Transformation: Y = sign(X) * log10(1 + |X| / C)
%
%   SYMLOG(AX, DIMS, C) transforms the specified axes ('x', 'y', or 'z') of
%   the graphics objects in AX to a symmetric logarithmic scale. It also
%   updates the axes ticks and labels to reflect the original values (e.g.,
%   10^1, -10^2) while maintaining the linear scaling around zero.
%
%   INPUTS:
%       (Mode 1: Data Transformation)
%       X       - (numeric) Input array to transform.
%       C       - (numeric, optional) Linear threshold parameter (default = 1).
%                 Values inside [-C, C] are compressed linearly.
%
%       (Mode 2: Axis Scaling)
%       AX      - (Axes Handle, optional) Target axes. Default: gca.
%       DIMS    - (char/string) Dimensions to transform: 'x', 'y', 'z', 'xy'...
%                 Default: 'y'.
%       C       - (numeric, optional) Linear threshold (default = 1).
%
%   OUTPUTS:
%       Y       - (numeric) Transformed data (only in Mode 1).
%
%   EXAMPLES:
%       % 1. Transform Data
%       y = symlog(x);
%
%       % 2. Transform Axis (Visual)
%       scatter(x, y);
%       symlog(gca, 'x');  % Apply value transformation + tick update
%
%   ALGORITHM:
%       The "Modified Log" transformation is used:
%           f(x) = log10(1 + |x|/C) * sign(x)
%       This function behaves like log10(|x|) for |x| >> C, and linearly
%       approximates x for |x| << C. It is continuous and smooth at 0.
%
%   See also: LOG10, SEMILOGX, SEMILOGY

%% ========================================================================
%  PARSE INPUTS
%  ========================================================================

% Check for Data Transformation Mode (Numeric Input)
% If first input is numeric and not a handle, or if it's a handle array but
% clearly data (unlikely for handle), treat as data.
if nargin >= 1 && isnumeric(varargin{1}) && ~all(ishghandle(varargin{1}(:)))
    % MODE 1: Data Transform
    x = varargin{1};
    if nargin >= 2
        C = varargin{2};
    else
        C = 1;
    end
    
    % Perform Transform
    varargout{1} = do_transform(x, C);
    return;
end

% MODE 2: Axis Scaling
args = varargin;
ax = [];
dims = '';
C = 1;

% 1. Extract Axis Handle
if ~isempty(args) && all(ishghandle(args{1}(:))) && strcmpi(get(args{1}(1), 'Type'), 'axes')
    ax = args{1};
    args(1) = [];
else
    ax = gca;
end

% 2. Extract Dimensions
if ~isempty(args) && (ischar(args{1}) || isstring(args{1}))
    dims = char(args{1});
    args(1) = [];
else
    dims = 'y';
end

% 3. Extract Constant C
if ~isempty(args) && isnumeric(args{1})
    C = args{1};
end

%% ========================================================================
%  AXIS TRANSFORMATION
%  ========================================================================

% Ensure Dimensions are lower case
dims = lower(dims);
validDims = {'x', 'y', 'z'};

for iDim = 1:length(validDims)
    d = validDims{iDim};
    
    % Check if this dimension is requested
    if contains(dims, d)
        apply_symlog_axis(ax, d, C);
    end
end


end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function y = do_transform(x, C)
    % Core mathematical transformation
    y = sign(x) .* log10(1 + abs(x) / C);
end

function x = do_inverse(y, C)
    % Inverse transformation: x = sign(y) * C * (10^|y| - 1)
    x = sign(y) .* C .* (10.^abs(y) - 1);
end

function apply_symlog_axis(ax, dim, C)
    % Apply transformation to a specific dimension of an axes
    
    % 1. Transform all graphic objects in the axes
    % We search for objects that have XData/YData/ZData properties.
    propName = [upper(dim) 'Data'];
    
    % Find all children (recursive)
    objs = findall(ax);
    
    % Filter for objects having the property
    hasProp = isprop(objs, propName);
    objs = objs(hasProp);
    
    for iObj = 1:length(objs)
        % Get raw data
        raw = get(objs(iObj), propName);
        
        % Avoid double-transforming? 
        % ideally check userdata or appdata, but simple approach:
        % Assume user calls this once per plotting command.
        
        if isnumeric(raw)
            trans = do_transform(raw, C);
            set(objs(iObj), propName, trans);
        end
        
        % Also handle BaseLine if it exists (for correlations, bars)
        if isprop(objs(iObj), 'BaseValue') && strcmpi(dim, 'y') 
            % Only typically relevant for Y (Bar plots)
            bv = get(objs(iObj), 'BaseValue');
            set(objs(iObj), 'BaseValue', do_transform(bv, C));
        end
    end
    
    % 2. Update Ticks to look nice
    % First, get the current limits in transformed space to know range
    % (We need to update limits if they were manual? No, MATLAB updates lims auto usually)
    % Force update of limits
    % axis(ax, 'auto'); % Can disturb other dims
    
    % Get new limits
    limName = [upper(dim) 'Lim'];
    lims = get(ax, limName);
    
    % Generate "Nice" ticks in LOG domain
    % Strategy: Powers of 10 in original domain -> transform to log domain
    
    % Inverse transform limits to find range in Real Units
    realLims = do_inverse(lims, C);
    
    % Determine powers of 10 to cover range
    % e.g. -100 to 100
    % Powers: ..., -100, -10, -1, 0, 1, 10, 100, ...
    
    maxPow = ceil(log10(max(abs(realLims)) + eps));
    % Generate powers: 0, 10^0, 10^1...
    pows = 0:maxPow;
    valsPos = 10.^pows;
    vals = unique([-fliplr(valsPos), 0, valsPos]);
    
    % Filter vals to be within/near limits (add some margin)
    % Actually, we can just show all within reason.
    
    % Transform these 'nice' values
    tickPos = do_transform(vals, C);
    
    % Filter ticks within visible range (optional, or let MATLAB clip)
    % Keep ticks that are within limits (relaxed)
    inRange = tickPos >= lims(1) & tickPos <= lims(2);
    
    % If too few ticks, maybe add intermediate?
    % For strict Matplotlib replication: Powers of 10.
    
    % Apply Ticks
    tickPosFinal = tickPos;
    
    % Generate Labels
    % 0 -> '0'
    % 10^k -> '10^k'
    % -10^k -> '-10^k'
    
    labels = cell(size(tickPosFinal));
    for i = 1:length(vals)
        v = vals(i);
        if v == 0
            labels{i} = '0';
        elseif abs(v) < 10 && abs(v) >= 1 && mod(v,1)==0
             % Integers < 10 (e.g. 1) -> Just '1' or '-1'? 
             % Matplotlib often uses 10^0. Let's use 10^x notation mostly.
             % User requested "tidy". '-10^1' is clear. '1' is clearer than '10^0'.
             labels{i} = sprintf('%g', v);
        elseif abs(log10(abs(v))) == floor(abs(log10(abs(v)))) 
            % Exact power of 10
            p = log10(abs(v));
            if v > 0
                labels{i} = sprintf('10^{%d}', p);
            else
                labels{i} = sprintf('-10^{%d}', p);
            end
        else
            % Just number
            labels{i} = sprintf('%g', v);
        end
    end
    
    % Set Ticks
    set(ax, [upper(dim) 'Tick'], tickPosFinal);
    set(ax, [upper(dim) 'TickLabel'], labels);
    
    % Add 'Mode: linear' indicator? No, keep clean.
end
