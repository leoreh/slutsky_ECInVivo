function plot_hypnogram(varargin)
% PLOT_HYPNOGRAM Plot a hypnogram as colored horizontal lines.
%
%   PLOT_HYPNOGRAM(varargin)
%
%   SUMMARY:
%       Draws one horizontal line per vigilance state across its bouts. Bouts
%       can be passed directly (boutTimes) or decoded from a per-bin labels
%       vector. All drawing targets the supplied axis (hAx); a figure is
%       created only when no axis is given.
%
%   INPUTS (Name-Value or positional):
%       labels      - (Vec)  1 x n integer state labels. See as_classify.
%       boutTimes   - (Cell) {nstates x 1} of [start end] matrices (units of
%                            the desired x-axis, e.g. seconds or hours).
%       sstates     - (Vec)  Selected states to mark (order sets overlap priority).
%       clr         - (Cell) Per-state RGB colors; if empty loads cfg.colors.
%       yshift      - (Num)  Scalar shift of line location ('lines' only). {1}
%       lWidth      - (Num)  Line width ('lines' only). {20}
%       hAx         - (Handle) Axis to draw into. If empty, a new figure is made.
%       style       - (Char) 'lines' {default} draws colored horizontal lines;
%                            'strip' draws a compact filled color strip (one
%                            patch per bout, tight y-band [0 1]).
%
%   DEPENDENCIES:
%       as_loadConfig, as_bouts (only when decoding labels).
%
%   HISTORY:
%       18 Jun 22 LH
%       22 Jun 26 axis-safe (all ops target hAx); house-style header

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addOptional(p, 'labels', [], @isnumeric);
addOptional(p, 'boutTimes', [], @iscell);
addOptional(p, 'sstates', [], @isnumeric);
addOptional(p, 'clr', []);
addOptional(p, 'yshift', 1, @isnumeric);
addOptional(p, 'lWidth', 20, @isnumeric);
addOptional(p, 'hAx', []);
addParameter(p, 'style', 'lines', @(x) any(strcmpi(char(x), {'lines', 'strip'})));

parse(p, varargin{:})
labels          = p.Results.labels;
boutTimes       = p.Results.boutTimes;
sstates         = p.Results.sstates;
clr             = p.Results.clr;
yshift          = p.Results.yshift;
lWidth          = p.Results.lWidth;
hAx             = p.Results.hAx;
style           = char(p.Results.style);

if isempty(boutTimes) && isempty(labels)
    error('must input boutTimes or labels')
end

if isempty(hAx)
    fh = figure;
    hAx = subplot(1, 1, 1);
end

%% ========================================================================
%  PREPARE DATA
%  ========================================================================

% state params
cfg = as_loadConfig();

% selected states. the order of sstates determines which state is shown in
% case of overlap (from merging nearby bouts). e.g. setting state 4 (nrem)
% after state 2 (qw) shows a shared bin as nrem.
if isempty(sstates)
    if isempty(labels)
        sstates = 1 : cfg.nstates;
    else
        sstates = unique(labels(~isnan(labels)));
    end
end

% colors for states
if isempty(clr)
    clr = cfg.colors(sstates);
end

% re-calc state bouts from labels
if isempty(boutTimes)
    bouts = as_bouts('labels', labels, ...
        'minDur', 5, 'interDur', 3, 'graphics', false);
    boutTimes = bouts.times;
end

%% ========================================================================
%  PLOT
%  ========================================================================

hold(hAx, 'on')
switch lower(style)
    case 'lines'
        % one colored horizontal line per state at the top of the axis
        yLimit = ylim(hAx);
        for istate = 1 : length(sstates)
            sbouts = boutTimes{sstates(istate)};
            if ~isempty(sbouts)
                plot(hAx, sbouts', yLimit(2) * yshift * ones(size(sbouts))', ...
                    'color', clr{istate}, 'LineWidth', lWidth, ...
                    'HandleVisibility', 'off')
            end
        end
        yLimit = ylim(hAx);
        ylim(hAx, [yLimit(1), yLimit(2) * yshift])

    case 'strip'
        % compact filled color strip: one patch per bout spanning the full
        % y-band, tight (no whitespace). yshift / lWidth are ignored here.
        for istate = 1 : length(sstates)
            sbouts = boutTimes{sstates(istate)};
            if ~isempty(sbouts)
                x1 = sbouts(:, 1)';
                x2 = sbouts(:, 2)';
                X = [x1; x2; x2; x1];
                Y = repmat([0; 0; 1; 1], 1, numel(x1));
                patch(hAx, X, Y, clr{istate}, 'EdgeColor', 'none', ...
                    'HandleVisibility', 'off')
            end
        end
        ylim(hAx, [0, 1])
end
set(hAx, 'ytick', [])
set(hAx, 'YColor', 'none')

end

% EOF
