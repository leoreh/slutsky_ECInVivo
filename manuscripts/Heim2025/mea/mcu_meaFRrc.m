function mcu_meaFRrc(v, alt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEA Firing Rate Recovery Correlation Analysis
%
% Inputs:
%   v - cell array of data structures (if empty, loads data automatically)
%   alt - analysis option:
%       1: Plot correlation vs burstiness (B>=2 and B>=4)
%       2: Plot correlation vs dISI
%       3: Plot correlation per file with burstiness distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data if not provided
if isempty(v)
    grps = {'mea_bac'; 'mea_mcuko'};
    vars = {'mea', 'frr', 'st_brst'};
    basepaths = mcu_basepaths(grps{1}); % Only control group
    v = {basepaths2vars('basepath', basepaths, 'vars', vars)};
end

% Run selected analysis
switch alt
    case 1
        plot_corr_vs_burstiness(v);
    case 2
        plot_corr_vs_disi(v);
        % case 3
        %     plot_corr_per_file(v);
    otherwise
        error('Invalid alt option. Use 1, 2, or 3.');
end

end

function plot_corr_vs_burstiness(v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot correlation vs burstiness (B>=2 and B>=4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get data from control group
frr = catfields([v{1}(:).frr], 1);
uGood = frr.uGood;
spktimes = catfields([v{1}(:).mea], 2, false);
spktimes = spktimes.spktimes(uGood);

% Set parameters
mdlPrfx = 'mdlF';
rcvFld = 'spkDfct';
rcvFld = 'rcvWork';
rcvData = frr.(mdlPrfx).(rcvFld)(uGood);
isiThr = logspace(log10(0.005), 0, 20)';

% Calculate B>=2 correlation
[rVal2, ~, bspks2] = mea_frrc(spktimes, rcvData,...
    'minSpks', 2, 'isiThr', isiThr,...
    'mfr', [], 'flgPlot', false, 'flgBin', false);

% Calculate B>=4 correlation
[rVal4, ~, bspks4] = mea_frrc(spktimes, rcvData,...
    'minSpks', 4, 'isiThr', logspace(log10(0.005), 0, 10)',...
    'mfr', [], 'flgPlot', false, 'flgBin', false);

% Plot B>=2 correlation
[hFig, hAx] = plot_frrc(rVal2, bspks2, [], isiThr, false);
hPlt = findobj(hAx, 'Type', 'line', '-and', 'LineStyle', '-');
hPlt(2).DisplayName = 'B ≥ 2 Spikes';
hPlt(1).Color = [0, 0.5, 0.5];
hAx.YAxis(2).Color = [0, 0.4, 0.4];
hAx.XAxis.Label.String = 'Burst Freq. (Hz)';

% Plot B>=4 correlation
xVal = 1 ./ logspace(log10(0.005), 0, 10)';
yyaxis(hAx, 'left');
hCorr = plot(xVal, rVal4, '--', 'Color', [0.7, 0.7, 0.7, 0.7], 'LineWidth', 2);
axis tight;
hCorr.DisplayName = 'B ≥ 4 Spikes';
hCorr.Color = repmat(0.5, 1, 3);

% Plot B>=4 BSpks
yyaxis right
ylim([0 1])
hBurst = plot(xVal, mean(bspks4, 1, 'omitnan'), '--', 'LineWidth', 2);
hBurst.LineStyle = '--';
hBurst.LineWidth = 2;
hBurst.HandleVisibility = 'off';
hBurst.Color = [0.1, 0.65, 0.65];

% Formatting
hLgd = legend([hPlt(2), hCorr], 'Location', 'best', 'Interpreter', 'tex');
hLgd.Box = 'on';

if strcmp(rcvFld, 'rcvWork')
    txtTtl = 'Recovery Ratio';
    hLgd.Location = 'northeast';
    yyaxis left
    ylim([-0.1 0.3]);
else
    txtTtl = 'Spike Deficit';
    hLgd.Location = 'northeast';
    yyaxis left
    ylim([-0.3 0.11]);
end
set(hAx, 'FontSize', 16);
hAx.YAxis(1).Label.FontSize = 20;
hAx.YAxis(2).Label.FontSize = 20;
hAx.Title.FontSize = 20;

title(hAx, txtTtl)
plot_axSize('hFig', hFig, 'szOnly', true,...
    'axShape', 'square', 'axHeight', 300);

% Save
fname = sprintf('MEA ~ Corr(B,%s)', rcvFld);
lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat'});

end

function plot_corr_vs_disi(v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot correlation vs dISI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get data
frr = catfields([v{1}(:).frr], 1);
uGood = frr.uGood;
spktimes = catfields([v{1}(:).mea], 2, false);
spktimes = spktimes.spktimes(uGood);

% Set parameters
mdlType = 'mdlF';
isiThr = logspace(log10(0.005), 0, 20)';

% Calculate correlation
rcvFld = 'rcvWork';
rcvData = frr.(mdlType).(rcvFld)(uGood);
[rVal, ~, bspks] = mea_frrc(spktimes, rcvData,...
    'minSpks', 2, 'isiThr', isiThr,...
    'mfr', [], 'flgPlot', false, 'flgBin', true);

% Plot correlation
[hFig, hAx] = plot_frrc(rVal, bspks, [], isiThr, true);

% Graphics
hAx.XAxis.Label.String = 'Burst_{ΔISI} Freq. (Hz)';
if strcmp(rcvFld, 'rcvWork')
    txtTtl = 'Recovery Ratio';
else
    txtTtl = 'Spike Deficit';
end

title(hAx, txtTtl)
plot_axSize('hFig', hFig, 'szOnly', false,...
    'axShape', 'square', 'axHeight', 300);

% Save
fname = sprintf('MEA ~ Corr(BdISI,%s)', rcvFld);
lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat'});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [hFig, hAx] = plot_frrc(rVal, bspks, mfr, isiThr, flgBin)
% Plot correlation and burstiness analysis

% Create figure
[hFig, hAx] = plot_axSize('szOnly', false, 'axShape',...
    'square', 'axHeight', 300);
xVal = 1 ./ isiThr;

% Plot correlation coefficients
hCorr = plot(xVal, rVal, 'k', 'LineWidth', 2);
xticks([0.01, 0.1, 1, 10, 100])
set(hAx, 'XScale', 'log')
set(hAx, 'XDir', 'reverse')
hAx.XAxis.TickLabels = compose('%g', hAx.XAxis.TickValues);
hold on
xlabel('Burst ISI Freq. (Hz)', 'Interpreter', 'tex')
ylabel('Correlation (ρ)', 'Interpreter', 'tex');

% Add mean firing rate line
if ~isempty(mfr)
    hMfr = xline(mfr, '--r', 'DisplayName', 'MFR', 'LineWidth', 1);
end
% Add horizontal line at correlation = 0
yline(0, '--k', 'LineWidth', 1, 'Alpha', 0.5, 'HandleVisibility', 'off');

if isempty(bspks)
    return
end

% Normalize
binWidths = diff([0; isiThr]);
if flgBin
    dataMat = bspks ./ binWidths';
    dataMat = dataMat ./ sum(dataMat, 2);
    yTxt = 'P(S\inB)_{ΔISI}';
    yLbl = 'P(S\inB) Prob. Density';
else
    dataMat = bspks;
    yTxt = 'P(S\inB)';
    yLbl = 'P(S\inB) Cum. Prob.';
end
dataMat(isinf(dataMat) | isnan(dataMat)) = 0;

% Add right y-axis for histogram data
yyaxis right
clr = [0, 0.4, 0.4];
hAx.YAxis(2).Color = clr;
ylabel(yLbl, 'Interpreter', 'tex', 'Rotation', 270)

% Plot burstiness using plot_stdShade
hBurst = plot_stdShade('dataMat', dataMat,...
    'xVal', xVal, 'hAx', hAx, ...
    'clr', clr, 'alpha', 0.5);
hBurst.LineStyle = '-';
hBurst.LineWidth = 2;
hBurst.DisplayName = yTxt;

% Format figure
plot_axSize('hFig', hFig, 'szOnly', false,...
    'axShape', 'square', 'axHeight', 300);

end












% function plot_corr_per_file(v)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot correlation per file with burstiness distribution
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % Set parameters
% mdlType = 'mdlF';
% rcvFld = 'spkDfct';
% flgBin = false;
% nFiles = length(v{1});
%
% % Get data
% frr = catfields([v{1}(:).frr], 1);
% uGood = frr.uGood;
% spktimesAll = catfields([v{1}(:).mea], 2, false).spktimes(uGood);
% rcvDataAll = frr.(mdlType).(rcvFld)(uGood);
% frBslAll = frr.(mdlType).frBsl(uGood);
%
% % Initialize arrays
% rVal = nan(20, nFiles);
% mfr = nan(1, nFiles);
% bspks = cell(1, nFiles);
%
% % Calculate per file
% for iFile = 1 : nFiles
%     frr = v{1}(iFile).frr;
%     uGood = frr.uGood;
%     spktimes = v{1}(iFile).mea.spktimes(uGood);
%     rcvData = frr.(mdlType).(rcvFld)(uGood);
%     mfr(iFile) = mean(frr.(mdlType).frBsl(uGood), 1, 'omitnan');
%
%     [rVal(:, iFile), ~, bspks{iFile}] = mea_frrc(spktimes, rcvData,...
%         'mfr', mfr(iFile), 'flgPlot', false, 'flgBin', flgBin);
% end
%
% % Define ISI thresholds
% isiThr = logspace(-2.3, 1, 20)';
%
% % Prepare burstiness data
% dataMat = zeros(size(rVal));
% for iFile = 1:nFiles
%     if ~isempty(bspks{iFile})
%         dataMat(:, iFile) = mean(bspks{iFile}, 1, 'omitnan');
%     end
% end
%
% % Plot using helper function
% [hFig, hAx] = plot_frrc(rVal, dataMat, mfr, isiThr, flgBin);

% end