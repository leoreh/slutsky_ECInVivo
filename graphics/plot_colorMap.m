function [hIm, hBar] = plot_colorMap(data, varargin)
% PLOT_COLORMAP Displays a 2D matrix as a scaled color image.
%
%   [hIm, hBar] = PLOT_COLORMAP(data, varargin)
%
%   SUMMARY:
%       In-repo replacement for FMAToolbox's PlotColorMap. Wraps IMAGESC
%       with defaults suited to neural heatmaps (PETHs, firing maps):
%       explicit x/y axis vectors, color limits that autoscale robustly
%       over NaN/Inf, 'normal' y-direction, out-facing ticks, and an
%       optional colorbar. NaN entries render transparent (axes background).
%
%   INPUTS:
%       data        - (mat)  [M x N] Matrix to display (rows -> y, cols -> x).
%                            N-D slices (e.g. [1 x E x B]) are squeezed to 2D.
%       varargin    - (param/value) Optional parameters:
%           'hAx'      - (axes)     Target axes {gca}.
%           'x'        - (vec)      N abscissae mapped to columns {1:N}.
%           'y'        - (vec)      M ordinates mapped to rows {1:M}.
%           'cutoffs'  - (vec)      [lo hi] color limits {[] -> autoscale}.
%           'clrMap'   - (char/mat) Colormap name or [K x 3] matrix {'turbo'}.
%           'flgBar'   - (log)      Draw a colorbar? {false}.
%           'ydir'     - (char)     'normal' {default} | 'reverse'.
%
%   OUTPUTS:
%       hIm         - (handle) Image object (from IMAGESC).
%       hBar        - (handle) Colorbar object, or [] when 'flgBar' is false.
%
%   See also: IMAGESC, PLOT_SPEC.
%
%   HISTORY:
%       Jul 2026 - Vendored from FMAToolbox PlotColorMap (M. Zugaro), rewritten.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'data', @isnumeric);
addParameter(p, 'hAx', [], @(x) isempty(x) || isgraphics(x));
addParameter(p, 'x', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'y', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'cutoffs', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
addParameter(p, 'clrMap', 'turbo', @(x) ischar(x) || isstring(x) || (isnumeric(x) && size(x, 2) == 3));
addParameter(p, 'flgBar', false, @islogical);
addParameter(p, 'ydir', 'normal', @(x) any(validatestring(x, {'normal', 'reverse'})));

parse(p, data, varargin{:});
hAx     = p.Results.hAx;
xVal    = p.Results.x;
yVal    = p.Results.y;
cutoffs = p.Results.cutoffs;
clrMap  = p.Results.clrMap;
flgBar  = p.Results.flgBar;
ydir    = p.Results.ydir;

if isempty(hAx), hAx = gca; end
if isstring(clrMap), clrMap = char(clrMap); end

%% ========================================================================
%  PREPARE DATA
%  ========================================================================

% Collapse singleton dims for N-D slices (e.g. [1 x nEvents x nBins]).
% Only squeeze true N-D arrays so a genuine [1 x N] row is left intact.
if ~ismatrix(data)
    data = squeeze(data);
end
if ~ismatrix(data)
    error('plot_colorMap:notMatrix', 'data must reduce to a 2-D matrix.');
end
[nRows, nCols] = size(data);

% Axis vectors: default to indices; warn (do not error) on length mismatch
if isempty(xVal)
    xVal = 1:nCols;
elseif numel(xVal) ~= nCols
    warning('plot_colorMap:xMismatch', ...
        'numel(x)=%d ~= %d columns; using column indices.', numel(xVal), nCols);
    xVal = 1:nCols;
end
if isempty(yVal)
    yVal = 1:nRows;
elseif numel(yVal) ~= nRows
    warning('plot_colorMap:yMismatch', ...
        'numel(y)=%d ~= %d rows; using row indices.', numel(yVal), nRows);
    yVal = 1:nRows;
end
xVal = xVal(:);
yVal = yVal(:);

% Color limits: explicit cutoffs else autoscale over finite entries
if isempty(cutoffs)
    finiteVals = data(isfinite(data));
    lo = min(finiteVals);
    hi = max(finiteVals);
else
    lo = cutoffs(1);
    hi = cutoffs(2);
end

% Guard degenerate ranges (empty, non-finite, or flat)
if isempty(lo) || isempty(hi) || ~isfinite(lo) || ~isfinite(hi)
    lo = 0; hi = 1;
elseif lo == hi
    hi = lo + 1;
end

%% ========================================================================
%  PLOT
%  ========================================================================
hIm = imagesc(hAx, xVal, yVal, data);

% Manual CLim (version-robust vs clim/caxis) and NaN transparency
set(hAx, 'CLim', [lo hi]);
set(hIm, 'AlphaData', ~isnan(data));

% Axes cosmetics
set(hAx, 'YDir', ydir, 'TickDir', 'out', 'Box', 'off');
axis(hAx, 'tight');
colormap(hAx, clrMap);

% Colorbar
if flgBar
    hBar = colorbar(hAx);
    set(hBar, 'TickDir', 'out', 'Box', 'off');
else
    hBar = [];
end

end     % EOF
