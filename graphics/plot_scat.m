function [hSct, hFit] = plot_scat(tbl, x, y, varargin)
% PLOT_SCAT Flexible scatter plot with grouping, regression, and stats.
%
%   PLOT_SCAT(TBL, X, Y, ...) creates a scatter plot of variable Y vs X.
%   INPUTS:
%       tbl         - (table) or empty [].
%       x, y        - (char) Column name in tbl OR (numeric) Vector.
%
%   OPTIONAL KEY-VALUE PAIRS:
%       'g'         - (char/vector) Grouping variable.
%       'c'         - (char/vector/matrix) Color.
%                     If Variable/Vector (N x 1): Maps values to color.
%                     If Matrix (k x 3): Maps groups to colors.
%                     If Single Color (char/1x3): Fixed color.
%       'sz'        - (char/scalar/vector) Marker size (Default: 20).
%       'alpha'     - (scalar) Marker transparency (0-1).
%       'marker'    - (char) Marker type (Default: 'o').
%       'fitType'   - (char) 'Linear', 'Ortho', or 'None' (Default).
%       'flgStats'  - (logical) Show correlation stats (Default: false).
%       'dispName'  - (char) Legend display name.
%       'hAx'       - (handle) Target axis (Default: gca).
%
%   OUTPUTS:
%       hSct, hFit
%
%   See also: PLOT_LINEREG

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @(x) istable(x) || isempty(x));
addRequired(p, 'x', @(x) ischar(x) || isstring(x) || isnumeric(x));
addRequired(p, 'y', @(x) ischar(x) || isstring(x) || isnumeric(x));

addParameter(p, 'g', [], @(x) ischar(x) || isstring(x) || isnumeric(x) || iscategorical(x));
addParameter(p, 'c', [], @(x) ischar(x) || isstring(x) || isnumeric(x));
addParameter(p, 'sz', 20, @(x) ischar(x) || isstring(x) || isnumeric(x));
addParameter(p, 'alpha', 0.6, @isnumeric);
addParameter(p, 'marker', 'o', @(x) ischar(x) || isstring(x));
addParameter(p, 'fitType', 'None', @ischar);
addParameter(p, 'flgStats', false, @islogical);
addParameter(p, 'dispName', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'hAx', [], @(x) isempty(x) || isgraphics(x));

parse(p, tbl, x, y, varargin{:});
args = p.Results;

if isempty(args.hAx), args.hAx = gca; end

%% ========================================================================
%  DATA PREPARATION
%  ========================================================================

% --- Unified Extraction ---
xData  = extractData(tbl, args.x);
yData  = extractData(tbl, args.y);
gData  = extractData(tbl, args.g);
szData = extractData(tbl, args.sz);
cData  = extractData(tbl, args.c);

% --- Grouping ---
if isempty(gData)
    gData = ones(size(xData)); 
    isGrp = false;
    uGrps = 1;
else
    isGrp = true;
    if islogical(gData) || isnumeric(gData), gData = categorical(gData); end
    uGrps = unique(gData);
end
nGrps = length(uGrps);

% --- Color Logic ---
% Modes: 'data' (per-point), 'map' (per-group from map), 'fixed' (single), 'auto' (lines)
cMode = 'auto'; 
grpColors = [];

if ~isempty(cData)
    if (isnumeric(cData) || iscategorical(cData) || isstring(cData)) && ...
            length(cData) == length(xData)
        cMode = 'data';
    elseif isnumeric(cData) && size(cData, 2) == 3
        % Colormap (Nx3 or 1x3)
        cMode = 'map';
        if size(cData, 1) < nGrps
             % Cycle if map is smaller than groups
             if size(cData, 1) == 1
                 grpColors = repmat(cData, nGrps, 1);
             else
                 grpColors = repmat(cData, ceil(nGrps/size(cData,1)), 1);
             end
        else
             grpColors = cData;
        end
    elseif ischar(cData) || isstring(cData)
        % Fixed Color Name (e.g. 'r', 'red')
        cMode = 'fixed';
        if isstring(cData) && isscalar(cData), cData = char(cData); end
        grpColors = repmat(cData, nGrps, 1);
    end
end

if strcmp(cMode, 'auto')
    grpColors = lines(nGrps);
end

%% ========================================================================
%  PLOTTING
%  ========================================================================

wasHold = ishold(args.hAx);
hold(args.hAx, 'on');

hSct = gobjects(nGrps, 1);
hFit = gobjects(nGrps, 1);

for iGrp = 1:nGrps
    curGrp = uGrps(iGrp);
    
    % Index
    if isGrp
        grpIdx = (gData == curGrp);
        grpName = string(curGrp);
    else
        grpIdx = true(size(xData));
        grpName = 'All';
    end
    
    xG = xData(grpIdx);
    yG = yData(grpIdx);
    
    if isempty(xG), continue; end
       
    % Size
    if isscalar(szData)
        sG = szData;
    else
        sG = szData(grpIdx);
    end
    
    % Prepare Scatter Arguments
    if strcmp(cMode, 'data')
        % scatter(x,y,sz,c, ...)
        cG = cData(grpIdx);
        scatterArgs = {sG, cG, 'filled'}; 
        fitClr = [0 0 0]; % Fit line black if points are multi-colored
    else
        % 'map', 'fixed', 'auto'
        % Use MarkerFaceColor
        cG = grpColors(iGrp, :);
        fitClr = cG;
        scatterArgs = {sG, 'MarkerFaceColor', cG}; 
    end
    
    % Label
    if isempty(args.dispName)
        dispLbl = grpName;
    else
        if isGrp && nGrps > 1
            dispLbl = sprintf('%s (%s)', args.dispName, grpName);
        else
            dispLbl = args.dispName;
        end
    end
    
    % Draw Scatter
    hSct(iGrp) = scatter(args.hAx, xG, yG, scatterArgs{:}, ...
        'Marker', args.marker, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', args.alpha, ...
        'DisplayName', dispLbl);
        
    % Fit Line
    if ~strcmpi(args.fitType, 'None')
        [hF, ~] = plot_lineReg(xG, yG, 'hAx', args.hAx, 'type', args.fitType, ...
            'clr', fitClr, 'flgTxt', args.flgStats);
        
        hFit(iGrp) = hF;
    end
    
    % Stats
    if args.flgStats
        [rho, pval] = corr(xG, yG, 'Type', 'Spearman', 'Rows', 'complete');
        if pval < 0.0001, pStr = 'p<0.0001'; else, pStr = sprintf('p=%.4f', pval); end
        statsStr = sprintf('%c=%.2f, %s', char(961), rho, pStr);
        
        curLbl = get(hSct(iGrp), 'DisplayName');
        set(hSct(iGrp), 'DisplayName', sprintf('%s, %s', curLbl, statsStr));
    end
end

if ~wasHold
    hold(args.hAx, 'off');
end

end     % PLOT_SCAT


%% ========================================================================
%  HELPER: EXTRACTDATA
%  ========================================================================

function data = extractData(tbl, val)
% EXTRACTDATA Unified extraction.
    if isempty(val)
        data = [];
        return;
    end
    
    if istable(tbl) && (ischar(val) || isstring(val))
        if ismember(val, tbl.Properties.VariableNames)
            data = tbl.(val);
            return;
        end
    end
    
    data = val;
end
