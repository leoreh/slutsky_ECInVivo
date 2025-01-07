function psd_compare(psdCell, faxis, grps, clr, snames, txt_yax, varargin)
% Plot PSD Comparison for multiple states and groups, highlighting significant clusters.
%
% INPUT
%   psdCell     Cell array containing PSD data [state x freq x subject].
%   faxis       Frequency axis (1 x freq).
%   grps        Group names or identifiers.
%   clr         Colors for each state.
%   snames      State names for plot titles.
%   txt_yax     Y-axis label.
%
% VARARGIN
%   nperm       Number of permutations for cluster-based permutation {1000}.
%   alpha       Significance threshold {0.05}.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse inputs
p = inputParser;
addParameter(p, 'nperm', 1000, @isnumeric);
addParameter(p, 'alpha', 0.05, @isnumeric);
parse(p, varargin{:});

nperm = p.Results.nperm;
alpha = p.Results.alpha;

nstates = size(psdCell{1}, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot PSD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMatlabGraphics(true)
fh = figure;
set(fh, 'WindowState', 'maximized');
tlayout = [1, nstates];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
title(th, 'PSD Analysis', 'interpreter', 'none', 'FontSize', 20)
set(fh, 'DefaultAxesFontSize', 16);

% psd per state
for istate = 1 : nstates
    axh = nexttile(th, istate, [1, 1]); cla; hold on

    % Plot PSD for each group
    for igrp = 1 : length(grps)
        clr_state = (igrp == 1) * clr{istate} + (igrp == 2) * clr{istate} * 0.5;
        psdMat = squeeze(psdCell{igrp}(istate, :, :));
        plot_stdShade('dataMat', psdMat, 'xVal', faxis, 'axh', axh, 'clr', clr_state, 'alpha', 0.5)
        set(gca, 'yScale', 'log', 'XScale', 'log')
        xlabel('Frequency [Hz]')
        ylabel(txt_yax)

        ymax = max(mean(psdMat, 2));
        ymin = min(mean(psdMat, 2));
        yLimit = [10^floor(log10(ymin)), 10^ceil(log10(ymax))];
        set(gca, 'YLim', yLimit);

        psdGrp{igrp} = psdMat;
    end

    % Statistical comparison to find significant clusters
    stat = stat_cluPermutation(psdGrp, faxis);
    if isfield(stat, 'posclusters')
        if ~isempty(stat.posclusters)
            pos_bins = stat.freq(stat.posclusterslabelmat > 0 & stat.posclusters(1).prob < stat.cfg.alpha);
            plot(pos_bins, repmat(yLimit(2), 1, length(pos_bins)), 'r*', 'MarkerSize', 8);
        end
    end
    if isfield(stat, 'negclusters')
        if ~isempty(stat.negclusters)
            neg_bins = stat.freq(stat.negclusterslabelmat > 0 & stat.negclusters(1).prob < stat.cfg.alpha);
            plot(neg_bins, repmat(yLimit(2), 1, length(neg_bins)), 'b*', 'MarkerSize', 8);
        end
    end

    title(axh, snames{istate})
    legend({'WT', 'MCU-KO'}, 'Location', 'southwest')
end

end
