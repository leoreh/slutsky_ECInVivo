function [hHist, hKDE, hStat] = plot_hist(tbl, val, varargin)
% PLOT_HIST Flexible histogram with grouping, KDE, and stats.
%
%   PLOT_HIST(TBL, VAL, ...) creates a histogram of variable VAL from TBL.
%   PLOT_HIST([], VAL, ...) allows passing a vector directly for VAL.
%   INPUTS:
%       tbl         - (table) or empty [].
%       val         - (char) Column name in tbl OR (numeric) Vector.
%
%   OPTIONAL KEY-VALUE PAIRS:
%       'g'         - (char/vector) Grouping variable.
%       'c'         - (char/vector/matrix) Color.
%                     If Matrix (k x 3): Maps groups to colors.
%                     If Single Color (char/1x3): Fixed color (replicated).
%       'bins'      - (numeric) Number of bins OR Vector of edges.
%       'orient'    - (char) 'vertical' (default) or 'horizontal'.
%       'scale'     - (char) 'linear' (default) or 'log'.
%       'norm'      - (char) 'pdf' (default), 'probability', 'count'.
%       'flgKDE'    - (logical) Overlay Kernel Density Estimate (Default: false).
%       'flgStat'   - (logical) Overlay Mean line (Default: true).
%       'alpha'     - (scalar) Face alpha (Default: 0.3).
%       'hAx'       - (handle) Target axis (Default: gca).
%
%   OUTPUTS:
%       hHist, hKDE, hStat
%
%   See also: PLOT_SCAT, HISTOGRAM

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @(x) istable(x) || isempty(x));
addRequired(p, 'val', @(x) ischar(x) || isstring(x) || isnumeric(x));

addParameter(p, 'g', [], @(x) ischar(x) || isstring(x) || isnumeric(x) || iscategorical(x));
addParameter(p, 'c', [], @(x) ischar(x) || isstring(x) || isnumeric(x));
addParameter(p, 'bins', [], @isnumeric);
addParameter(p, 'orient', 'vertical', @(x) any(validatestring(x, {'vertical', 'horizontal'})));
addParameter(p, 'scale', 'linear', @(x) any(validatestring(x, {'linear', 'log'})));
addParameter(p, 'norm', 'pdf', @ischar);
addParameter(p, 'flgKDE', false, @islogical);
addParameter(p, 'flgStat', true, @islogical);
addParameter(p, 'alpha', 0.3, @isnumeric);
addParameter(p, 'hAx', [], @(x) isempty(x) || isgraphics(x));

parse(p, tbl, val, varargin{:});
args = p.Results;

if isempty(args.hAx), args.hAx = gca; end

%% ========================================================================
%  DATA PREPARATION
%  ========================================================================

% --- Unified Extraction ---
vData = extractData(tbl, args.val);
gData = extractData(tbl, args.g);

% --- Grouping ---
if isempty(gData)
    gData = ones(size(vData)); 
    isGrp = false;
    uGrps = 1;
else
    isGrp = true;
    if islogical(gData) || isnumeric(gData), gData = categorical(gData); end
    uGrps = unique(gData);
end
nGrps = length(uGrps);

% --- Color Logic ---
cIn = args.c;
grpColors = [];

if isempty(cIn)
    % Auto (lines)
    grpColors = lines(nGrps);
else
    % User Provided Color(s)
    if ischar(cIn) || isstring(cIn)
        % Fixed Color Name -> Convert to RGB
        try 
            cVal = validatecolor(cIn);
            % Replicate for all groups
            grpColors = repmat(cVal, nGrps, 1);
        catch
            % Fallback
            grpColors = lines(nGrps);
        end
    elseif isnumeric(cIn)
        % Matrix or Vector
        if size(cIn, 2) == 3
            % It is a colormap (N x 3) or single RGB color (1 x 3)
            if size(cIn, 1) == 1
                % Single color -> Replicate
                grpColors = repmat(cIn, nGrps, 1);
            elseif size(cIn, 1) < nGrps
                % Cycle map if smaller than groups
                grpColors = repmat(cIn, ceil(nGrps/size(cIn,1)), 1);
            else
                % Use map as is (or truncated)
                grpColors = cIn;
            end
        else
            % Invalid numeric format for histogram -> Fallback
            grpColors = lines(nGrps);
        end
    end
end

% --- Binning ---
bins = args.bins;
if isempty(bins)
    validV = vData(~isnan(vData) & ~isinf(vData));
    
    if strcmp(args.scale, 'log')
        validV = validV(validV > 0);
        mn=min(validV);
        mx=max(validV);
        bins = logspace(log10(mn), log10(mx), 30);
    else
        mn=min(validV);
        mx=max(validV);
        bins = linspace(mn, mx, 30);
    end
end

%% ========================================================================
%  PLOTTING
%  ========================================================================

wasHold = ishold(args.hAx);
hold(args.hAx, 'on');

hHist = gobjects(nGrps, 1);
hKDE  = gobjects(nGrps, 1);
hStat = gobjects(nGrps, 1);

for iGrp = 1:nGrps
    curGrp = uGrps(iGrp);
    
    % Index
    if isGrp
        grpIdx = (gData == curGrp);
        grpName = string(curGrp);
    else
        grpIdx = true(size(vData));
        grpName = 'All';
    end
    
    subV = vData(grpIdx);
    
    % Filter Valid
    subV = subV(~isnan(subV) & ~isinf(subV));
    if strcmp(args.scale, 'log')
        subV = subV(subV > 0);
    end
    
    if isempty(subV), continue; end
    
    % Style
    col = grpColors(iGrp, :);
    
    % Histogram
    hHist(iGrp) = histogram(args.hAx, subV, ...
        'BinEdges', bins, ...
        'FaceColor', col, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', args.alpha, ...
        'Normalization', args.norm, ...
        'Orientation', args.orient, ...
        'DisplayName', grpName);
    
    % KDE
    if args.flgKDE && length(subV) > 1
        hKDE(iGrp) = plotKDE(args.hAx, subV, bins, col, args.orient, args.scale, args.norm);
        set(hKDE(iGrp), 'HandleVisibility', 'off');
    end
    
    % Stats (Mean)
    if args.flgStat
        statVal = mean(subV);
        statArgs = {'LineStyle', '-', 'Color', col, 'LineWidth', 2, 'HandleVisibility', 'off'};
        
        if strcmp(args.orient, 'vertical')
            yl = ylim(args.hAx);
            hStat(iGrp) = plot(args.hAx, [statVal statVal], yl, statArgs{:});
        else
            xl = xlim(args.hAx);
            hStat(iGrp) = plot(args.hAx, xl, [statVal statVal], statArgs{:});
        end
    end
end

if ~wasHold
    hold(args.hAx, 'off');
end

end     % PLOT_HIST


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
        
        if strcmp(normType, 'probability')
            % Probability Adaptation
            binW = mean(diff(log10(edges))); 
            f = f_log * binW;
        else
             % PDF
             f = f_log; 
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
