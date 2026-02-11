function [hHist, hKDE, hMean] = plot_hist(tbl, val, varargin)
% PLOT_HIST Flexible histogram with grouping, KDE, and stats.
%
%   PLOT_HIST(TBL, VAL, ...) creates a histogram of variable VAL from TBL.
%
%   PLOT_HIST([], VAL, ...) allows passing a vector directly for VAL.
%
%   INPUTS:
%       tbl         - (table) or empty.
%       val         - (char) Column name in tbl OR (numeric) Vector.
%
%   OPTIONAL KEY-VALUE PAIRS:
%       'g'         - (char) Grouping variable name in tbl OR (vector) grouping.
%       'c'         - (char) Color variable name OR (matrix) RGB colors.
%       'bins'      - (numeric) Number of bins OR Vector of edges.
%       'orientation' - (char) 'vertical' (default) or 'horizontal'.
%       'scale'     - (char) 'linear' (default) or 'log'.
%       'normalization' - (char) 'pdf' (default), 'probability', 'count'.
%       'showKDE'   - (logical) Overlay Kernel Density Estimate (Default: false).
%       'showMean'  - (logical) Overlay Mean line (Default: false).
%       'showMedian'- (logical) Overlay Median line (Default: false).
%       'alpha'     - (scalar) Face alpha (Default: 0.3).
%       'hAx'       - (handle) Target axis (Default: gca).
%
%   OUTPUTS:
%       hHist       - (handle) Histogram object(s).
%       hKDE        - (handle) KDE line object(s).
%       hMean       - (handle) Mean/Median line object(s).
%
%   See also: PLOT_SCAT, HISTOGRAM

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @(x) istable(x) || isempty(x));
addRequired(p, 'val', @(x) ischar(x) || isstring(x) || isnumeric(x));

% Optional
addParameter(p, 'g', [], @(x) ischar(x) || isstring(x) || isnumeric(x) || iscategorical(x));
addParameter(p, 'c', [], @(x) ischar(x) || isstring(x) || isnumeric(x));
addParameter(p, 'bins', [], @isnumeric);
addParameter(p, 'orientation', 'vertical', @(x) any(validatestring(x, {'vertical', 'horizontal'})));
addParameter(p, 'scale', 'linear', @(x) any(validatestring(x, {'linear', 'log'})));
addParameter(p, 'normalization', 'pdf', @ischar);
addParameter(p, 'showKDE', false, @islogical);
addParameter(p, 'showMean', false, @islogical);
addParameter(p, 'showMedian', false, @islogical);
addParameter(p, 'alpha', 0.3, @isnumeric);
addParameter(p, 'hAx', [], @(x) isempty(x) || isgraphics(x));

parse(p, tbl, val, varargin{:});

hAx = p.Results.hAx;
if isempty(hAx), hAx = gca; end

gIn = p.Results.g;
cIn = p.Results.c;
bins = p.Results.bins;
orient = p.Results.orientation;
scl = p.Results.scale;
normType = p.Results.normalization;
flgKDE = p.Results.showKDE;
flgMean = p.Results.showMean;
flgMed = p.Results.showMedian;
alp = p.Results.alpha;

%% ========================================================================
%  DATA PREPARATION
%  ========================================================================

if istable(tbl)
    if ischar(val) || isstring(val), vData = tbl.(val); else, vData = val; end
    
    if isempty(gIn)
        gData = ones(height(tbl), 1);
        uGrps = 1;
        isGrp = false;
    else
        if ischar(gIn) || isstring(gIn)
            gData = tbl.(gIn);
        else
            gData = gIn;
        end
        if islogical(gData), gData = categorical(gData); end
        if ~iscategorical(gData) && ~isnumeric(gData), gData = categorical(gData); end
        uGrps = unique(gData);
        isGrp = true;
    end
    
else
    vData = val;
    if isempty(gIn)
        gData = ones(length(vData), 1);
        uGrps = 1;
        isGrp = false;
    else
        gData = gIn;
        if islogical(gData), gData = categorical(gData); end
        if ~iscategorical(gData) && ~isnumeric(gData), gData = categorical(gData); end
        uGrps = unique(gData);
        isGrp = true;
    end
end

% Colors
cMap = [];
if isnumeric(cIn) && size(cIn, 2) == 3
    cMap = cIn;
elseif (ischar(cIn) || isstring(cIn))
    try
        cMap = validatecolor(cIn);
    catch
        % ignore if column name? logic here similar to plot_scat
    end
end

%% ========================================================================
%  PLOTTING
%  ========================================================================

hold(hAx, 'on');

hHist = gobjects(length(uGrps), 1);
hKDE  = gobjects(length(uGrps), 1);
hMean = gobjects(length(uGrps), 1);

% Color cycling
if isempty(cMap)
    defColors = lines(length(uGrps));
else
    if size(cMap, 1) == 1 && length(uGrps) > 1
        defColors = repmat(cMap, length(uGrps), 1);
    elseif size(cMap, 1) < length(uGrps)
        defColors = repmat(cMap, ceil(length(uGrps)/size(cMap,1)), 1);
    else
        defColors = cMap;
    end
end

% Binning Strategy (Global if not provided)
if isempty(bins)
    % Calculate global bins to ensure alignment
    validV = vData(~isnan(vData) & ~isinf(vData));
    if strcmp(scl, 'log')
        validV = validV(validV > 0);
        if isempty(validV), mn=0.1; mx=1; else, mn=min(validV); mx=max(validV); end
        % Logspace bins
        bins = logspace(log10(mn), log10(mx), 30);
    else
        if isempty(validV), mn=0; mx=1; else, mn=min(validV); mx=max(validV); end
        bins = linspace(mn, mx, 30);
    end
end

% Log Scale Normalization Adjustment
% If log scale, 'pdf' is often misleading visually unless using 'probability'.
% But users might request parameter. We will respect input but adjust KDE scaling.
% Logic from tblGUI_scatHist (plot_histPdf)
if strcmp(scl, 'log') && strcmp(normType, 'pdf')
    % Often we want 'probability' for log visual
    % But let's stick to user request. If they say pdf, we give pdf.
end


for iG = 1:length(uGrps)
    currGrp = uGrps(iG);
    idx = (gData == currGrp);
    subV = vData(idx);
    
    % Filter valid
    subV = subV(~isnan(subV) & ~isinf(subV));
    if strcmp(scl, 'log')
        subV = subV(subV > 0); 
    end
    
    if isempty(subV), continue; end
    
    col = defColors(iG, :);
    
    % Name
    if isGrp
        if iscategorical(currGrp) || isstring(currGrp), nm = string(currGrp); else, nm = num2str(currGrp); end
    else
        nm = 'All';
    end
    
    % Histogram
    hHist(iG) = histogram(hAx, subV, ...
        'BinEdges', bins, ...
        'FaceColor', col, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', alp, ...
        'Normalization', normType, ...
        'Orientation', orient, ...
        'DisplayName', nm);
    
    % KDE
    if flgKDE && length(subV) > 1
        hKDE(iG) = plotKDE(hAx, subV, bins, col, orient, scl, normType);
        set(hKDE(iG), 'HandleVisibility', 'off');
    end
    
    % Mean / Median
    if flgMean || flgMed
        if flgMean, statVal = mean(subV); lSty = '-'; end
        if flgMed,  statVal = median(subV); lSty = '--'; end
        
        % Plot Line
        if strcmp(orient, 'vertical')
            yl = ylim(hAx);
            hMean(iG) = plot(hAx, [statVal statVal], yl, 'LineStyle', lSty, ...
                'Color', col, 'LineWidth', 2, 'HandleVisibility', 'off');
        else
            xl = xlim(hAx);
            hMean(iG) = plot(hAx, xl, [statVal statVal], 'LineStyle', lSty, ...
                'Color', col, 'LineWidth', 2, 'HandleVisibility', 'off');
        end
    end
    
end

hold(hAx, 'off');

end

%% ========================================================================
%  HELPER: KDE
%  ========================================================================
function hLine = plotKDE(ax, val, edges, c, orient, scale, normType)
    
    % Determine points
    if strcmp(scale, 'log')
        logMin = log10(min(edges));
        logMax = log10(max(edges));
        pts    = logspace(logMin, logMax, 200);
        
        % KDE on log data
        try
             [f_log, ~] = ksdensity(log10(val), log10(pts));
        catch
             hLine = gobjects(0); return; 
        end
        
        % Scale transformation logic
        % If we want 'probability' mass match:
        % The histogram height is Prob/BinWidth (if pdf) or Prob (if probability).
        % Standard ksdensity returns PDF (integral = 1).
        
        if strcmp(normType, 'probability')
            % To match probability bars (which are integral over bin), we multiply by average bin width
            % But on log scale, bin width varies. This is tricky.
            % Adaptation from tblGUI_scatHist strategy:
            binW = mean(diff(log10(edges))); 
            % f_log is dP/d(log10x). 
            % We want approx Prob per bin ~ dP * BinWidth
            f = f_log * binW;
        else
             % Assume 'pdf'. But PDF on log axis is different.
             % Basic visualization: just plot the log-density shape scaled to fit
             % For now, stick to the robust tblGUI method if probability, else raw PDF
             f = f_log; 
             % This might need refinement for 'pdf' on log axis to match histogram height perfectly
             % but is visually consistent with density.
        end
    else
        pts = linspace(min(edges), max(edges), 200);
        try
            [f, ~] = ksdensity(val, pts);
        catch
            hLine = gobjects(0); return;
        end
        
        if strcmp(normType, 'probability')
             binW = mean(diff(edges));
             f = f * binW;
        end
    end
    
    if strcmp(orient, 'vertical')
        hLine = plot(ax, pts, f, 'Color', c, 'LineWidth', 2);
    else
        hLine = plot(ax, f, pts, 'Color', c, 'LineWidth', 2);
    end
end
