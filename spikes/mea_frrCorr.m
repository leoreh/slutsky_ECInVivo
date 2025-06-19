function [rVal, pVal, bspks] = mea_frrCorr(spktimes, rcvData, varargin)
% MEAFRRCORR Calculates correlations between burstiness and recovery metrics.
%
% SUMMARY:
% This function analyzes the relationship between neuronal burstiness and
% firing rate recovery parameters. It systematically varies burst detection
% parameters (ISI threshold) to find optimal correlations with recovery metrics.
% The analysis can be run in two modes based on 'flgBin':
%   - Binned (flgBin=true): Burstiness is calculated for ISIs within discrete bins.
%   - Cumulative (flgBin=false): Burstiness is calculated for all ISIs up to a given threshold.
%
% The analysis includes:
%   1. Burst detection using variable ISI thresholds/bins
%   2. Calculation of burstiness metrics (spike fraction in bursts)
%   3. Correlation analysis with recovery metrics (spike deficit, recovery gain, etc.)
%   4. Generation of a summary plot showing correlation patterns and burstiness distribution.
%
% INPUT (Required):
%   spktimes      - Cell array. spktimes{i} contains spike times (s) for neuron i (already filtered to uGood).
%   rcvData       - Vector of recovery metric values (already filtered to uGood).
%
% INPUT (Optional Key-Value Pairs):
%   winLim        - [start end] time window to analyze [s] {[0, 70*60]}.
%   isiThr        - Vector of ISI thresholds to test [s] {logspace(-2.3, 1, 20)'}.
%   flgPlot       - Logical flag to generate summary plot {true}.
%   flgBin        - Logical flag for binned (true) vs. cumulative (false) analysis {false}.
%   thrSpks       - Minimum number of spikes required per unit {0}.
%   minSpks       - Minimum spikes per burst {2}.
%   mfr           - Mean firing rate for plotting {[]}.
%
% OUTPUT:
%   rVal          - Correlation coefficients (nIsiThr x 1).
%   pVal          - P-values for correlations (nIsiThr x 1).
%   bspks         - Spike fraction in bursts (nUnits x nIsiThr).

% -------------------------------------------------------------------------
% INPUT PARSING
% -------------------------------------------------------------------------

p = inputParser;
addRequired(p, 'spktimes', @iscell);
addRequired(p, 'rcvData', @isnumeric);
addParameter(p, 'winLim', [0, 70 * 60], @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'isiThr', logspace(-2.3, 1, 20)', @isnumeric);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgBin', false, @islogical);
addParameter(p, 'thrSpks', 0, @(x) isnumeric(x) && isscalar(x) && x>=0);
addParameter(p, 'minSpks', 2, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'mfr', [], @isnumeric);

parse(p, spktimes, rcvData, varargin{:});
winLim = p.Results.winLim;
isiThr = p.Results.isiThr;
flgPlot = p.Results.flgPlot;
flgBin = p.Results.flgBin;
thrSpks = p.Results.thrSpks;
minSpks = p.Results.minSpks;
mfr = p.Results.mfr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA PREPARATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepares spike times and recovery metrics for correlation analysis:
% 1. Filter spike times to analysis window
% 2. Ensure minimum spike count requirements

% Filter spike times to analysis window
spktimes = cellfun(@(x) x(x >= winLim(1) & x <= winLim(2)),...
    spktimes, 'UniformOutput', false);

% Ensure minimum spikes
nSpks = cellfun(@length, spktimes);
spktimes = spktimes(nSpks >= thrSpks);
rcvData = rcvData(nSpks >= thrSpks);
nUnits = length(spktimes);

nThr = length(isiThr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BURSTINESS CALCULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates burstiness based on specified ISI thresholds/bins.
% The calculation method is determined by the `flgBin` parameter.

% Calculate burstiness for each ISI threshold/bin
bspks = zeros(nUnits, nThr);
isiThrEdges = [0; isiThr];
for iThr = 1:nThr
    for iUnit = 1:nUnits
        spks = spktimes{iUnit};
        if length(spks) < minSpks
            continue;
        end
        isi = diff(spks);

        % Detect bursts using binary2bouts
        if flgBin
            % ISIs within a specific bin: (previous thr, current thr]
            isBurstIsi = isi > isiThrEdges(iThr) & isi <= isiThrEdges(iThr+1);
        else
            % ISIs less than or equal to the threshold (cumulative)
            isBurstIsi = isi <= isiThr(iThr);
        end
        
        [idxBrsts, ~] = binary2bouts('vec', isBurstIsi, ...
            'minDur', minSpks-1, 'flgPrnt', false);
        
        if ~isempty(idxBrsts)
            % Calculate fraction spikes in bursts
            nBspks = diff(idxBrsts, 1, 2) + 1;
            bspks(iUnit, iThr) = sum(nBspks) / nSpks(iUnit);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRELATION ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates correlations between burstiness and recovery metrics:
% 1. Compute Spearman correlations with recovery metric for all ISI thresholds at once
% 2. Store results in correlation matrices

% Calculate correlation for all ISI thresholds at once using pairwise complete
[rVal, pVal] = corr(bspks, rcvData, 'type', 'Spearman', 'Rows', 'pairwise');

% Plotting
if flgPlot
    frrCorr_plot(bspks, rVal, isiThr, mfr, flgBin);
end

end     % EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: FRR_CORR_PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function frrCorr_plot(bspks, rVal, isiThr, mfr, flgBin)
% Create figure
[hFig, hAx] = plot_axSize('szOnly', false, 'axShape', 'square', 'axHeight', 300);
xVal = 1 ./ isiThr;

% Plot correlation coefficients
hCorr = plot(xVal, rVal, 'k', 'LineWidth', 2.5);
xticks([0.01, 0.1, 1, 10, 100])
set(hAx, 'XScale', 'log')
set(hAx, 'XDir', 'reverse')
hAx.XAxis.TickLabels = compose('%g', hAx.XAxis.TickValues);
hold on
xlabel('Burst ISI Freq. (Hz)', 'Interpreter', 'tex')
if flgBin
    yTxt = 'ρ(P(S\inB)_{ΔISI}, hDfct)';
    yTxt = 'Corr(SpkB, FR_{Deficit})';
else
    yTxt = 'Corr(SpkB, FR_{Deficit})';
end
hCorr.DisplayName = yTxt;
ylabel('Correlation (ρ)', 'Interpreter', 'tex');

% Add horizontal line at correlation = 0
yline(0, '--k', 'LineWidth', 1, 'Alpha', 0.5);

% ------------------------------------------------------------------------
% Burstiness PDF / CDF

% Normalize
binWidths = diff([0; isiThr]);
if flgBin
    dataMat = bspks ./ binWidths';
    dataMat = dataMat ./ sum(dataMat, 2);
    yTxt = 'P(S\inB)_{ΔISI}';
    yLbl = 'Probability Density';
else
    dataMat = bspks;
    yTxt = 'P(S\inB)';
    yLbl = 'Cum. Probability';
end
dataMat(isinf(dataMat) | isnan(dataMat)) = 0; 

% Add right y-axis for histogram data
yyaxis right
clr = [0, 0.4, 0.4];
hAx.YAxis(2).Color = clr;
ylabel(yLbl, 'Interpreter', 'tex', 'Rotation', 270)
% hAx.YAxis(2).TickLabels = {};
% hAx.YAxis(1).TickLabelFormat = '%.1f';

% Plot burstiness using plot_stdShade
hBurst = plot_stdShade('dataMat', dataMat,...
    'xVal', xVal, 'axh', hAx, ...
    'clr', clr, 'alpha', 0.5);
hBurst.LineStyle = '-';
hBurst.LineWidth = 1;
hBurst.DisplayName = yTxt;

% Add mean firing rate line
if ~isempty(mfr)
    hMfr = xline(mfr, '--r', 'DisplayName', 'MFR', 'LineWidth', 1);
end

% Create legend
% hLgd = legend([hCorr, hBurst], 'Location', 'best', 'Interpreter', 'tex');
% hLgd.Box = 'on';

% Format figure
plot_axSize('hFig', hFig, 'szOnly', false,...
    'axShape', 'square', 'axHeight', 300);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REMARKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For minSpks = 2, the "fraction of spikes in a burst" metric is a
% mathematical tautology (approximately linear transformation) of the ISI's
% cumulative distribution function evaluated at isiThr. It provides no new
% information about the spike train's structure that is not already present
% in the first-order ISI distribution itself.
% spkFrac ~ 2 * isiHist, where isiHist:
% isiHist = histcounts(isi, 'BinEdges', [0, isiThr'],...
%     'Normalization', 'cdf') * 2;
% The insistance on "spkFrac" is to maintain compatibility with different
% minSpks

% Since isiThr is logspaced, to plot the probability of a given spkFrac it
% must be divided by the bin width.