%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORR VS. BURSTINESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data 
iGrp = 1; 
grps = {'mea_bac'; 'mea_mcuko'};
vars = {'mea', 'frr', 'st_brst'};
basepaths = mcu_sessions(grps{iGrp});
v = basepaths2vars('basepath', basepaths, 'vars', vars);
nFiles = length(basepaths);

% Select param
winLim = [0, 70 * 60];     
mdlType = 'mdlF';
rcvFld = 'spkDfct';
% rcvFld = 'rcvGain';

% All files
clear rVal mfr brsty
frr = catfields([v(:).frr], 1);
uGood = frr.uGood;
spktimes = catfields([v(:).mea], 2, false);
spktimesAll = spktimes.spktimes(uGood);
rcvDataAll = frr.(mdlType).(rcvFld)(uGood);
frBslAll = frr.(mdlType).frBsl(uGood);
isiThr = logspace(log10(0.005), 0, 10)';
mfr = mean(frBslAll, 1, 'omitnan');

mea_frrCorr(spktimesAll, rcvDataAll,...
    'minSpks', 2, 'isiThr', logspace(log10(0.005), 0, 20)',...
    'mfr', [], 'flgPlot', true, 'flgBin', false);

hFig = gcf;
hAx = gca;
yyaxis left
hPlt = get(hAx, 'Children');
hPlt = hPlt(2);
hPlt.DisplayName = 'B ≥ 2 Spikes';

[rVal, ~, brsty] = mea_frrCorr(spktimesAll, rcvDataAll,...
    'minSpks', 4, 'isiThr', isiThr,...
    'mfr', [], 'flgPlot', false, 'flgBin', false);

hold on 
yyaxis right
% Plot burstiness using plot_stdShade
xVal = 1 ./ isiThr;
hBurst = plot_stdShade('dataMat', brsty,...
    'xVal', xVal, 'axh', hAx, ...
    'clr', [0, 0.7, 0.7], 'alpha', 0.5);
hBurst.LineStyle = '-';
hBurst.LineWidth = 1;
hBurst.HandleVisibility = 'off';
ylim([0 1])

yyaxis left
hCorr = plot(xVal, rVal, 'Color', [0.6, 0.6, 0.6, 0.6], 'LineWidth', 2.5);
ylim([-0.3 0.12])
% ylim([-0.05 0.22])
hCorr.DisplayName = 'B ≥ 4 Spikes';

% Create legend
hLgd = legend([hPlt, hCorr], 'Location', 'best', 'Interpreter', 'tex');
hLgd.Box = 'on';
hLgd.Location = 'northeast';

% Format figure
plot_axSize('hFig', hFig, 'szOnly', false,...
    'axShape', 'square', 'axWidth', 300);



fname = 'MEA ~ Corr(B,Gain)';
lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat'});




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dISI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select param
winLim = [0, 70 * 60];     
mdlType = 'mdlF';

% All files
frr = catfields([v(:).frr], 1);
uGood = frr.uGood;
spktimes = catfields([v(:).mea], 2, false);
spktimes = spktimes.spktimes(uGood);
frBsl = frr.(mdlType).frBsl(uGood);
isiThr = logspace(log10(0.005), 0, 10)';
mfr = mean(frBsl, 1, 'omitnan');

rcvFld = 'rcvGain';
rcvData = frr.(mdlType).(rcvFld)(uGood);
mea_frrCorr(spktimesAll, rcvData,...
    'minSpks', 2, 'isiThr', isiThr,...
    'mfr', mfr, 'flgPlot', true, 'flgBin', true);

hFig = gcf;
yyaxis left
ylim([-0.3, 0.4])
hAx = gca;
hPlt = get(hAx, 'Children');
hPlt = hPlt(2);
hPlt.DisplayName = 'ρ(B_{ΔISI}, Gain)';

rcvFld = 'spkDfct';
rcvData = frr.(mdlType).(rcvFld)(uGood);
[rVal, ~, ~] = mea_frrCorr(spktimesAll, rcvDataAll,...
    'minSpks', 2, 'isiThr', isiThr,...
    'mfr', [], 'flgPlot', false, 'flgBin', true);

hCorr = plot(xVal, rVal, 'Color', 'k', 'LineWidth', 2.5);
hCorr.LineStyle = '-';
hCorr.Color = [0.6, 0.6, 0.6, 0.6];
hCorr.DisplayName = 'ρ(B_{ΔISI}, Deficit)';

yyaxis right
hAx = gca;
hMfr = get(hAx, 'Children');
hMfr = hMfr(1);
hMfr.DisplayName = 'MFR';
ylim([0, 0.35])
yticks([0 : 0.1 : 0.5])

% Create legend
hLgd = legend([hPlt, hCorr, hMfr], 'Location', 'best', 'Interpreter', 'tex');
hLgd.Box = 'on';
hLgd.Location = 'northeast';

% Format figure
plot_axSize('hFig', hFig, 'szOnly', false,...
    'axShape', 'square', 'axWidth', 300);

fname = 'MEA ~ Corr(BdISI)';
lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PER FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flgRand = false;
flgBin = false;
nRand = 60;         % Number of units to randomly sample

% Per file analysis
rVal = nan(20, nFiles);
mfr = nan(1, nFiles);
brsty = cell(1, nFiles);

for iFile = 1 : nFiles
    if flgRand     
        uRand = randsample(sum(uGood), nRand);
        spktimes = spktimesAll(uRand);
        rcvData = rcvDataAll(uRand);
        mfr(iFile) = mean(frBslAll(uRand), 1, 'omitnan');
    else
        frr = v(iFile).frr;
        uGood = frr.uGood;
        spktimes = v(iFile).mea.spktimes(uGood);
        rcvData = frr.(mdlType).(rcvFld)(uGood);
        mfr(iFile) = mean(frr.(mdlType).frBsl(uGood), 1, 'omitnan');
    end

    [rVal(:, iFile), ~, brsty{iFile}] = mea_frrCorr(spktimes, rcvData,...
        'mfr', mfr(iFile), 'flgPlot', false, 'flgBin', flgBin);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create figure
[hFig, hAx] = plot_axSize('szOnly', false, 'axShape', 'square', 'axHeight', 300);

% Define ISI thresholds (same as in mea_frrCorr)
isiThr = logspace(-2.3, 1, 20)';
xVal = 1 ./ isiThr;

% Plot correlation coefficients using plot_stdShade
plot_stdShade('dataMat', rVal, 'xVal', xVal, 'axh', hAx, ...
    'clr', [0, 0, 0], 'alpha', 0.3);

% Format x-axis
xticks([0.01, 0.1, 1, 10, 100]);
set(hAx, 'XScale', 'log');
set(hAx, 'XDir', 'reverse');
xlabel('Burst Frequency (Hz)');
ylabel('ρ(Brsty(cum), Recovery)');

% Add horizontal line at correlation = 0
yline(0, '--k', 'LineWidth', 1, 'Alpha', 0.5);

% ------------------------------------------------------------------------
% Burstiness PDF / CDF

% Add right y-axis for histogram data
yyaxis right;
clr = [0, 0.5, 0.5];
hAx.YAxis(2).Color = clr;
ylabel('Brst Prob. Dist.');

% Prepare burstiness data for plotting
binWidths = diff([0; isiThr]);

% Normalize burstiness data across files
dataMat = zeros(size(rVal));
for iFile = 1:nFiles
    if ~isempty(brsty{iFile})
        if flgBin
            tempData = brsty{iFile} ./ binWidths';
        else
            tempData = brsty{iFile};
        end
        tempData = tempData ./ sum(tempData, 2);
        tempData(isinf(tempData) | isnan(tempData)) = 0;
        dataMat(:, iFile) = mean(tempData, 1, 'omitnan');
    end
end

% Plot burstiness using plot_stdShade
plot_stdShade('dataMat', dataMat, 'xVal', xVal, 'axh', hAx, ...
    'clr', clr, 'alpha', 0.5);

% Add vertical dashed red line for each MFR
xline(mfr, '--r', 'LineWidth', 1, 'Alpha', 0.5);

% Format figure
plot_axSize('hFig', hFig, 'szOnly', false,...
    'axShape', 'square', 'axHeight', 300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fname = lme_frml2char('MEA ~ Corr(B,S)', 'rmRnd', true);
% lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat'});










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRELATION BSPKS VS FR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data 
iGrp = 1; 
grps = {'mea_bac'; 'mea_mcuko'};
vars = {'mea', 'frr', 'st_brst'};
basepaths = mcu_sessions(grps{iGrp});
v = basepaths2vars('basepath', basepaths, 'vars', vars);
nFiles = length(basepaths);

% Select param
winLim = [0, 70 * 60];     
mdlType = 'mdlF';
rcvFld = 'spkDfct';
% rcvFld = 'rcvGain';

% All files
frr = catfields([v(:).frr], 1);
uGood = frr.uGood;

spktimes = catfields([v(:).mea], 2, false);
spktimes = spktimes.spktimes(uGood);
frBsl = frr.(mdlType).frBsl(uGood);
isiThr = logspace(log10(0.005), 0, 20)';
mfr = mean(frBsl, 1, 'omitnan');

rcvData = frBsl;

mea_frrCorr(spktimes, rcvData,...
    'minSpks', 2, 'isiThr', isiThr,...
    'mfr', [], 'flgPlot', true, 'flgBin', false);
