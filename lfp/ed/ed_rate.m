function edRate = ed_rate(ed, varargin)
% ED_RATE Epileptiform-discharge rate over time.
%
%   edRate = ED_RATE(ed, varargin)
%
%   SUMMARY:
%       Bins accepted discharge peak times into a rate time series via
%       times2rate. Optionally splits counts by vigilance state and plots.
%
%   INPUTS:
%       ed          - (Struct) Requires .peakTime; uses .accepted / .state if present.
%       varargin    - Parameter/Value pairs:
%           'binsize' - (Num)  Bin width [s]. {60}
%           'winCalc' - (Vec)  [M x 2] windows to bin within. {full extent}
%           'flgPlot' - (Log)  Plot the rate trace? {false}
%           'basepath'- (Char) For the plot title. {pwd}
%
%   OUTPUTS:
%       edRate      - (Struct):
%           .rate        (1 x nBins) rate [Hz]
%           .tstamps     (1 x nBins) bin centers [s]
%           .binEdges    (nBins x 2) bin edges [s]
%           .binsize     scalar [s]
%           .stateNames  (cell) categories present (if ed.state given)
%           .stateCounts (vec)  accepted-event count per state (if ed.state given)
%
%   DEPENDENCIES:
%       times2rate.
%
%   HISTORY:
%       Created: 22 Jun 2026

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'ed', @isstruct);
addParameter(p, 'binsize', 60, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'winCalc', [], @isnumeric);
addParameter(p, 'flgPlot', false, @islogical);
addParameter(p, 'basepath', pwd, @ischar);

parse(p, ed, varargin{:});
ed       = p.Results.ed;
binsize  = p.Results.binsize;
winCalc  = p.Results.winCalc;
flgPlot  = p.Results.flgPlot;
basepath = p.Results.basepath;

%% ========================================================================
%  RATE
%  ========================================================================

peakTime = ed.peakTime(:);
acc = true(size(peakTime));
if isfield(ed, 'accepted') && numel(ed.accepted) == numel(peakTime)
    acc = logical(ed.accepted(:));
end

[r, binEdges, binCents] = times2rate(peakTime(acc), ...
    'binsize', binsize, 'winCalc', winCalc, 'c2r', true);

edRate = struct();
edRate.rate     = r(:)';
edRate.tstamps  = binCents(:)';
edRate.binEdges = binEdges;
edRate.binsize  = binsize;

% Optional per-state counts
if isfield(ed, 'state') && ~isempty(ed.state)
    s = ed.state(acc);
    edRate.stateNames  = categories(s);
    edRate.stateCounts = countcats(s);
end

%% ========================================================================
%  PLOT
%  ========================================================================

if flgPlot
    [~, basename] = fileparts(basepath);
    fh = figure('Name', [basename, ' - ED rate'], 'NumberTitle', 'off');
    plot(edRate.tstamps / 3600, edRate.rate, 'k', 'LineWidth', 1);
    xlabel('Time (h)'); ylabel('ED rate (Hz)');
    title([basename, ' - ED rate'], 'Interpreter', 'none');
    axis tight; box off;
end

end     % EOF
