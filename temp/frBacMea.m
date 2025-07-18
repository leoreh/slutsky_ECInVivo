clear params
params.winLim = [0 60 * 60];        % Analysis window [s]
params.binSize = 0.002;             % 1ms bins
params.gkHw = 0.012;                % 12ms sigma
params.winStpr = 1.0;               % 1s window
params.nShuffles = 1000;             % Number of shuffles
params.flg_par = false;
params.spkLim = 2000;

basepaths = [mcu_sessions('wt_bsl'), mcu_sessions('mcu_bsl')];
nFiles = length(basepaths);
vars = {'spikes'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

tic
for iFile = 1 : nFiles

    basepath = basepaths{iFile};
    cd(basepath);
    [prc] = prCoupling(v(iFile).spikes.times, params, 'flgSave', true);
    prCoupling_plot(prc, 'basepath', basepath, 'flgSaveFig', true);

end
toc






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Organize files
grps = {'wt_bsl'; 'mcu_bsl'};

clear grppaths
for igrp = 1 : length(grps)
    grppaths{igrp} = string(mcu_sessions(grps{igrp})');
end

% -------------------------------------------------------------------------
% FR per unit, WT vs MCU for RS vs FS 
frml = 'PRC ~ Group * UnitType + (1|Mouse)';
varFld = 'z0';

% get data
[lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varFld', varFld);

% run lme
contrasts = 'all';
contrasts = [1 : 6]; 
[lmeStats, lmeCfg] = lme_analyse(lmeData, lmeCfg, 'contrasts', contrasts);

% plot
hFig = lme_plot(lmeData, lmeCfg.lmeMdl, 'ptype', 'bar', 'axShape', 'square'); 

% Update labels
hAx = gca;
ylabel(hAx, 'Firing Rate (Hz)')
xlabel(hAx, '')
title(hAx, '')
hAx.Legend.Location = 'northwest';

% add significance lines
barIdx = {[1, 2], [3, 4]};
barLbl = {'NS', '****'};
plot_sigLines(hAx, barIdx, barLbl, 'lineGap', 0.15)
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'square', 'axHeight', 300)

% save
fname = lme_frml2char(frml, 'rmRnd', true,...
    'sfx', ['_altClassify', num2str(altClassify)]);
lme_save('fh', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeData', lmeData, 'lmeStats', lmeStats, 'lmeCfg', lmeCfg)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRELATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grab data
frml = 'PRC ~ Group * UnitType + (1|Mouse)';
varFld = 'z0';
[lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varFld', varFld);
frml = 'Burst ~ Group * UnitType + (1|Mouse)';
varFld = 'royer';
[brstData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varFld', varFld);
frml = 'FR ~ Group * UnitType + (1|Mouse)';
[frData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false);

% Find shared UnitIDs
sharedIDs = intersect(brstData.UnitID, lmeData.UnitID);
[~, sharedIdx] = ismember(sharedIDs, brstData.UnitID);

% Select only shared units
%

% Create aligned tables with only shared units
prcTableShared = prcTable(idxPRC, :);
burstTableShared = burstTable(idxBurst, :);

% Verify alignment
assert(all(prcTableShared.UnitID == burstTableShared.UnitID), 'Tables not properly aligned');

% Plot baseline vs recovery firing rate
figure;

grpStr = {'Control', 'MCU-KO'};
unitStr = {'pPYR', 'pINT'};

iGrp = 1;
grpIdx = lmeData.Group == grpStr{iGrp};
for iUnit = 1 : 2

    unitIdx = lmeData.UnitType == unitStr{iUnit};
    plotIdx = unitIdx & grpIdx;
    x = prc.z0;

    subplot(2,2,iUnit);
plot_scatterCorr(prc.stpr0, frr.bslFr, ...
    'xLbl', 'Baseline Firing Rate (Hz)', 'yLbl', 'Recovery Firing Rate (Hz)');

end


















iFile = 1;
spktimes = v(iFile).mea.spktimes;
spktimes = mea.spktimes;

profile on;
[prc] = prCoupling(v(iFile).spikes.times(1 : 10), params, 'flgSave', true);
profile viewer;
profile off;

[prc] = prCoupling(v(iFile).spikes.times, params, 'flgSave', true);


frr = mea_frRecovery(spktimes, 'winLim', [0, 9 * 60 * 60]);


% Plot baseline vs recovery firing rate
figure;
subplot(1,2,1);
plot_scatterCorr(prc.stpr0, frr.bslFr, ...
    'xLbl', 'Baseline Firing Rate (Hz)', 'yLbl', 'Recovery Firing Rate (Hz)');

subplot(1,2,2);
xLbl = '';
yLbl = '';
plot_scatterCorr(prc.stpr0, frr.frRecov, 'flgOtl', true, ...
    'xLbl', xLbl, 'yLbl', yLbl);






















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flg_load = false;

if flg_load

    % load all files
    vars = {'mea', 'fr'};
    basepaths = mcu_sessions('mea_bac');
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);
    nFiles = length(basepaths);

    % Limit spktimes to first hour
    win = [0, 70 * 60];

    popC = popFRcoupling(spikes.times, 'win', win);

    for iFile = 1 : nFiles
        spktimes = v(iFile).mea.spktimes;

        fr = calc_fr(spktimes, 'basepath', pwd,...
            'graphics', false, 'binsize', 30, 'saveVar', false, 'forceA', true,...
            'smet', 'GK', 'winBL', win, 'winCalc', [0, Inf]);

        figure
        subplot(2,1,1)
        plot(fr.tstamps / 60 / 60, fr.strd)
        subplot(2,1,2)
        plot(fr.tstamps / 60 / 60, fr.norm)
    end

    spktimes = {};
    for iFile = 1 : nFiles
        spktimes = [spktimes, v(iFile).mea.spktimes];
    end

    fr = calc_fr(spktimes, 'basepath', pwd,...
        'graphics', false, 'binsize', 30, 'saveVar', false, 'forceA', true,...
        'smet', 'GK', 'winBL', win, 'winCalc', [0, Inf]);

    % spike timing metrics
    st = spktimes_metrics('spktimes', spktimes, 'sunits', [],...
        'bins', win, 'flg_force', true, 'flg_save', false, 'flg_all', false);

    % brst (mea)
    brst = spktimes_meaBrst(spktimes, 'binsize', [], 'isiThr', 0.005,...
        'minSpks', 2, 'flgSave', false, 'flgForce', true, 'bins', win);

    % Calculate fraction of ISIs below 10ms (100 Hz) for each cell
    bspks = zeros(1, length(spktimes));
    isiThr = 0.006;
    for iunit = 1 : length(spktimes)
        spku = spktimes{iunit}(spktimes{iunit} >= win(1) & spktimes{iunit} <= win(2));
        isi = diff(spku);
        bspks(iunit) = sum(isi < isiThr) / length(isi);
    end

end

% close all

% -------------------------------------------------------------------------
% inspect units
% figure;
% histogram(st.mizuseki)

% [~, otl] = rmoutliers(fr.fanoFactor);
% otl = fr.fanoFactor > 0.5;
% figure;
% plot(fr.strd(otl, :)')
% st.royer(otl)
% bspks(otl)

% params
brsty = brst.bspks;
flg_log = false;

% bad units
spkThr = 300;
nspks = cellfun(@(x) sum(InIntervals(x, win)), spktimes, 'uni', true);
badIdx = isnan(brsty) | isinf(brsty) | nspks < spkThr;

fprintf('nan: %d\n', sum(isnan(brsty)))
fprintf('zero: %d\n', sum(brsty == 0))
fprintf('inf: %d\n', sum(isinf(brsty)))
fprintf('nspks: %d\n', sum(nspks < spkThr))
fprintf('bad: %d\n', sum(badIdx))

% figure; plot(fr.strd(badIdx, :)')

brsty = brsty(~badIdx);
frt = fr.strd(~badIdx, :);

% Convert timestamps to hours
tstamps = fr.tstamps / 3600;  % Simplified conversion

% Define time intervals (in seconds, converted to hours in calculation)
tidx = [0, 1*3600; 2*3600, 3*3600; 8*3600, 9*3600] / 3600;

%% Calculate mean firing rates for intervals
% Initialize array for mean firing rates (units x intervals)
ufr = zeros(size(frt, 1), size(tidx, 1));

% Calculate mean firing rate for each interval
for iInterval = 1:size(tidx, 1)
    idx = tstamps >= tidx(iInterval,1) & tstamps <= tidx(iInterval,2);
    ufr(:,iInterval) = mean(frt(:,idx), 2);
end

%% Calculate normalized differences
% Raw differences
d12 = abs(ufr(:,2) - ufr(:,1));  % Absolute difference 2-1
d23 = abs(ufr(:,3) - ufr(:,2));  % Absolute difference 3-2
d13 = abs(ufr(:,3) - ufr(:,1));  % Absolute difference 3-1 (total change)

% Normalized absolute differences as percent of baseline
nd12 = 100 * (d12 ./ ufr(:,1));     % Normalized absolute change 2-1
nd23 = 100 * (d23 ./ ufr(:,2));     % Normalized absolute change 3-2
nd13 = 100 * (d13 ./ ufr(:,1));     % Normalized absolute change 3-1 (total change)

d12 = log10(d12);
d23 = log10(d23);
d13 = log10(d13);
nd12 = log10(nd12);
nd23 = log10(nd23);
nd13 = log10(nd13);

%% Calculate deviations from population mean
% Calculate mean FR across population for each interval
mfr = mean(ufr, 1);  % Mean across units for each interval

% Calculate absolute deviation from population mean as percent of mean
ufrLog = log10(ufr + eps);
mfr = mean(ufrLog, 1);
fr_dev = 100 * (abs((ufrLog - mfr) ./ mfr));  % Absolute percent deviation from mean

% figure;
% histogram(nd12)

%% Plot relationships
figure('Position', [100 100 1200 500]);

% Plot burst rate vs FR deviation
x = brsty; y = fr_dev(:,3);
subplot(2,3,1)
mea_plotCorr(brsty, fr_dev(:,3), flg_log);
xlabel('Burstiness');
ylabel('|MFR Deviation| (% of mean)');
% set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

% fr
subplot(2,3,4)
mea_plotCorr(ufrLog(:,1), fr_dev(:,3), flg_log);
xlabel('BSL FR');
ylabel('|MFR Deviation| (% of mean)');
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')

% Plot burst rate vs FR deviation
subplot(2,3,2)
mea_plotCorr(brsty, nd13, flg_log);
xlabel('Burstiness');
ylabel('|FR Change| (%BSL)');
set(gca, 'XScale', 'log')

% fr
subplot(2,3,5)
mea_plotCorr(ufrLog(:,1), nd13, flg_log);
xlabel('BSL FR');
ylabel('|FR Change| (%BSL)');
% set(gca, 'XScale', 'log')

% fr
subplot(2,3,3)
mea_plotCorr(ufrLog(:,1), brsty, flg_log);
xlabel('BSL FR');
ylabel('Burstiness');
set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')

%% Analyze recovery to baseline for each unit
% Smooth firing rates using moving average
smoothFR = smoothdata(frt, 2, 'gaussian', 60);  % Smooth across time dimension with 30-bin window

% Define time windows (in hours)
baselineWin = [0, 1];  % 0-60 min for baseline
perturbWin = [1.5, 2];   % Perturbation window
recoveryWin = [2, 9];    % 2-9 hours for recovery period

% Find indices for different periods
baseIdx = tstamps >= baselineWin(1) & tstamps <= baselineWin(2);
perturbIdx = tstamps >= perturbWin(1) & tstamps <= perturbWin(2);
recovIdx = tstamps >= recoveryWin(1) & tstamps <= recoveryWin(2);

% Parameters for baseline recovery detection
windowSize = 20;  % Number of bins to check for stability
baselineThresh = 0.2;  % FR must be within 20% of baseline
minStableTime = 5 / 60;    % Must be stable for at least 5 min
recoveryThresh = 0.5;   % Units that don't recover above 50% of baseline are considered non-recovering

% Initialize arrays for results
nUnits = size(smoothFR, 1);
recoveryTimes = zeros(nUnits, 1);
baselineFR = zeros(nUnits, 1);
minFR = zeros(nUnits, 1);
recoveryType = zeros(nUnits, 1);  % 1: recovered to baseline, 2: overshot, 3: failed to recover
recoveryQuality = zeros(nUnits, 1);  % How close to baseline the unit recovers

% For each unit
for iUnit = 1:nUnits
    % Get baseline and recovery data
    baseFR = smoothFR(iUnit, baseIdx);
    recovFR = smoothFR(iUnit, recovIdx);
    recovTime = tstamps(recovIdx);

    % Calculate baseline and minimum FR
    baselineFR(iUnit) = mean(baseFR);
    minFR(iUnit) = min(recovFR);

    % Calculate final recovery quality
    finalFR = mean(recovFR(end-windowSize:end));
    recoveryQuality(iUnit) = finalFR / baselineFR(iUnit);

    % Classify recovery type
    if recoveryQuality(iUnit) < recoveryThresh
        recoveryType(iUnit) = 3;  % Failed to recover
    elseif recoveryQuality(iUnit) > 1.2  % More than 20% above baseline
        recoveryType(iUnit) = 2;  % Overshot
    else
        recoveryType(iUnit) = 1;  % Recovered to baseline
    end

    % Find when FR first reaches and stabilizes at baseline
    if recoveryType(iUnit) ~= 3  % Only look for recovery time if unit recovered
        stableFound = false;

        % First find all points that are at baseline level
        atBaseline = abs(recovFR - baselineFR(iUnit)) / baselineFR(iUnit) < baselineThresh;

        % Find continuous periods at baseline
        for i = 1:(length(recovFR) - windowSize)
            % Check if all points in window are at baseline
            if all(atBaseline(i:i+windowSize-1))
                % Check if this stability lasts for minStableTime
                if recovTime(i+windowSize-1) - recovTime(i) >= minStableTime
                    recoveryTimes(iUnit) = recovTime(i);
                    stableFound = true;
                    break;
                end
            end
        end

        if ~stableFound
            recoveryTimes(iUnit) = NaN;  % No stable period found
        end
    else
        recoveryTimes(iUnit) = NaN;  % No recovery time for non-recovering units
    end
end

% Plot results
figure('Position', [100 100 1200 800]);

% Plot distribution of recovery times (only for recovering units)
subplot(2,2,1)
histogram(recoveryTimes(~isnan(recoveryTimes)), 20)
xlabel('Time to Baseline (hours)')
ylabel('Number of Units')
title('Distribution of Recovery Times')

% Plot recovery quality
subplot(2,2,2)
hold on
histogram(recoveryQuality(recoveryType == 1), 20, 'FaceColor', 'b')
histogram(recoveryQuality(recoveryType == 2), 20, 'FaceColor', 'g')
histogram(recoveryQuality(recoveryType == 3), 20, 'FaceColor', 'r')
xlabel('Recovery Quality (Final/Baseline)')
ylabel('Number of Units')
title('Recovery Quality by Type')
legend('Recovered', 'Overshot', 'Failed to Recover')

% Plot example recovery curves
figure('Position', [100 100 800 400]);
hold on
% Plot 5 units with most typical recovery (closest to median recovery time)
medianRecovTime = median(recoveryTimes(~isnan(recoveryTimes)));
[~, idx] = sort(abs(recoveryTimes - medianRecovTime));
for i = 1:5
    unit = idx(i);
    if ~isnan(recoveryTimes(unit))
        plot(tstamps, smoothFR(unit,:), 'LineWidth', 1)
        xline(recoveryTimes(unit), '--', sprintf('%.1f h', recoveryTimes(unit)))
    end
end
xlabel('Time (hours)')
ylabel('Firing Rate (Hz)')
title('Example Recovery Curves')

%% Calculate recovery slope (3-4 hours)
% Find indices for 3-4 hour window
slopeWin = [2, 3];  % hours
slopeIdx = tstamps >= slopeWin(1) & tstamps <= slopeWin(2);

% Initialize array for slopes
recoverySlopes = zeros(nUnits, 1);

% Calculate slope for each unit
for iUnit = 1:nUnits
    % Get FR data for slope window
    slopeFR = smoothFR(iUnit, slopeIdx);
    slopeTime = tstamps(slopeIdx);

    % Fit line to get slope
    p = polyfit(slopeTime, slopeFR, 1);
    recoverySlopes(iUnit) = p(1);  % Slope in Hz/hour
end

%% Analyze relationships with burstiness
figure('Position', [100 100 1200 800]);

% Get valid indices for recovery times
validIdx = ~isnan(recoveryTimes);

% 1. Burstiness vs Recovery Time
subplot(2,3,1)
mea_plotCorr(log10(brsty(validIdx) + eps), recoveryTimes(validIdx), false);
xlabel('log10(Burstiness)')
ylabel('Recovery Time (hours)')
title('Burstiness vs Recovery Time')

% 2. Burstiness vs Recovery Quality
subplot(2,3,2)
mea_plotCorr(log10(brsty + eps), recoveryQuality, false);
xlabel('log10(Burstiness)')
ylabel('Recovery Quality (Final/Baseline)')
title('Burstiness vs Recovery Quality')
yline(1, '--', 'Baseline')

% 3. Burstiness vs Baseline FR
subplot(2,3,3)
mea_plotCorr(log10(brsty + eps), log10(baselineFR), false);
xlabel('log10(Burstiness)')
ylabel('log10(Baseline FR)')
title('Burstiness vs Baseline FR')

% 4. Burstiness vs Recovery Slope
subplot(2,3,4)
mea_plotCorr(log10(brsty + eps), recoverySlopes, false);
xlabel('log10(Burstiness)')
ylabel('Recovery Slope (Hz/hour)')
title('Burstiness vs Recovery Slope')

% 5. Recovery Slope vs Recovery Time
subplot(2,3,5)
mea_plotCorr(recoverySlopes(validIdx), recoveryTimes(validIdx), false);
xlabel('Recovery Slope (Hz/hour)')
ylabel('Recovery Time (hours)')
title('Recovery Slope vs Time')

% 6. Recovery Slope vs Recovery Quality
subplot(2,3,6)
mea_plotCorr(recoverySlopes, recoveryQuality, false);
xlabel('Recovery Slope (Hz/hour)')
ylabel('Recovery Quality (Final/Baseline)')
title('Recovery Slope vs Quality')
yline(1, '--', 'Baseline')

% Print statistical relationships
fprintf('\nRelationships with Burstiness:\n')
[r, p] = corrcoef(log10(brsty(validIdx) + eps), recoveryTimes(validIdx));
fprintf('Correlation with recovery time: r=%.2f, p=%.3f\n', r(1,2), p(1,2))

[r, p] = corrcoef(log10(brsty + eps), recoveryQuality);
fprintf('Correlation with recovery quality: r=%.2f, p=%.3f\n', r(1,2), p(1,2))

[r, p] = corrcoef(log10(brsty + eps), log10(baselineFR));
fprintf('Correlation with baseline FR: r=%.2f, p=%.3f\n', r(1,2), p(1,2))

fprintf('\nRelationships with Recovery Slope:\n')
[r, p] = corrcoef(log10(brsty + eps), recoverySlopes);
fprintf('Correlation with burstiness: r=%.2f, p=%.3f\n', r(1,2), p(1,2))

[r, p] = corrcoef(recoverySlopes(validIdx), recoveryTimes(validIdx));
fprintf('Correlation with recovery time: r=%.2f, p=%.3f\n', r(1,2), p(1,2))

[r, p] = corrcoef(recoverySlopes, recoveryQuality);
fprintf('Correlation with recovery quality: r=%.2f, p=%.3f\n', r(1,2), p(1,2))

% Plot example recovery curves with slope window highlighted
figure('Position', [100 100 800 400]);
hold on
% Plot 5 units with most typical recovery (closest to median recovery time)
medianRecovTime = median(recoveryTimes(~isnan(recoveryTimes)));
[~, idx] = sort(abs(recoveryTimes - medianRecovTime));
for i = 1:5
    unit = idx(i);
    if ~isnan(recoveryTimes(unit))
        plot(tstamps, smoothFR(unit,:), 'LineWidth', 1)
        xline(recoveryTimes(unit), '--', sprintf('%.1f h', recoveryTimes(unit)))
        % Plot slope window
        xline(slopeWin(1), ':', 'Slope window')
        xline(slopeWin(2), ':', '')
        % Plot fitted line for slope
        slopeFR = smoothFR(unit, slopeIdx);
        slopeTime = tstamps(slopeIdx);
        p = polyfit(slopeTime, slopeFR, 1);
        yfit = polyval(p, slopeTime);
        plot(slopeTime, yfit, 'r--', 'LineWidth', 1)
    end
end
xlabel('Time (hours)')
ylabel('Firing Rate (Hz)')
title('Example Recovery Curves with Slope Window')

%% Fit exponential recovery model
% Define model function: FR(t) = A + B*exp(-(t-t0)/τ) for t > t0
% where A = final steady state, B = drop magnitude, τ = time constant, t0 = drop time

% Initialize arrays for model parameters
nUnits = size(frt, 1);
modelParams = zeros(nUnits, 4);  % [A, B, τ, t0]
modelR2 = zeros(nUnits, 1);

% Define time windows
baseWin = [0, 1];  % Baseline window
dropWin = [1.5, 2];  % Drop window
recovWin = [2, 9];  % Recovery window

% Get indices for different periods
baseIdx = tstamps >= baseWin(1) & tstamps <= baseWin(2);
dropIdx = tstamps >= dropWin(1) & tstamps <= dropWin(2);
recovIdx = tstamps >= recovWin(1) & tstamps <= recovWin(2);

% For each unit
for iUnit = 1:nUnits
    % Get baseline and recovery data
    baseFR = frt(iUnit, baseIdx);
    recovFR = frt(iUnit, recovIdx);
    recovTime = tstamps(recovIdx);

    % Initial parameter estimates
    A0 = mean(baseFR);  % Final steady state ≈ baseline
    B0 = A0 - min(recovFR);  % Drop magnitude
    t0 = dropWin(1);  % Drop time
    tau0 = 1;  % Initial time constant guess

    % Define model function
    modelFun = @(params, t) params(1) + params(2)*exp(-(t-params(4))/params(3));

    % Define objective function (sum of squared errors)
    objFun = @(params) sum((recovFR - modelFun(params, recovTime)).^2);

    % Set parameter bounds
    lb = [0, 0, 0.1, dropWin(1)];  % Lower bounds
    ub = [inf, inf, 10, dropWin(2)];  % Upper bounds

    % Fit model using fmincon
    options = optimoptions('fmincon', 'Display', 'off');
    params = fmincon(objFun, [A0, B0, tau0, t0], [], [], [], [], lb, ub, [], options);

    % Store parameters
    modelParams(iUnit,:) = params;

    % Calculate R-squared
    yfit = modelFun(params, recovTime);
    SSresid = sum((recovFR - yfit).^2);
    SStotal = (length(recovFR)-1) * var(recovFR);
    modelR2(iUnit) = 1 - SSresid/SStotal;
end

% Plot example fits
figure('Position', [100 100 1200 400]);

% Plot 5 units with best fits
[~, bestIdx] = sort(modelR2, 'descend');
for i = 1:5
    unit = bestIdx(i);
    subplot(1,5,i)
    hold on

    % Plot raw data
    plot(tstamps, frt(unit,:), 'k.', 'MarkerSize', 4)

    % Plot model fit
    recovTime = tstamps(recovIdx);
    yfit = modelParams(unit,1) + modelParams(unit,2)*exp(-(recovTime-modelParams(unit,4))/modelParams(unit,3));
    plot(recovTime, yfit, 'r-', 'LineWidth', 2)

    % Plot baseline
    yline(modelParams(unit,1), '--', sprintf('A=%.1f', modelParams(unit,1)))

    % Add parameter text
    text(0.05, 0.95, sprintf('τ=%.2f\nB=%.1f\nR²=%.2f', ...
        modelParams(unit,3), modelParams(unit,2), modelR2(unit)), ...
        'Units', 'normalized', 'VerticalAlignment', 'top')

    xlabel('Time (hours)')
    ylabel('Firing Rate (Hz)')
    title(sprintf('Unit %d', unit))
end

% Analyze relationships with model parameters
figure('Position', [100 100 1200 400]);

% 1. Burstiness vs Time Constant
subplot(1,3,1)
mea_plotCorr(log10(brsty + eps), modelParams(:,3), false);
xlabel('log10(Burstiness)')
ylabel('Recovery Time Constant (hours)')
title('Burstiness vs τ')

% 2. Burstiness vs Drop Magnitude
subplot(1,3,2)
mea_plotCorr(log10(brsty + eps), modelParams(:,2), false);
xlabel('log10(Burstiness)')
ylabel('Drop Magnitude (Hz)')
title('Burstiness vs B')

% 3. Time Constant vs Drop Magnitude
subplot(1,3,3)
mea_plotCorr(modelParams(:,3), modelParams(:,2), false);
xlabel('Time Constant (hours)')
ylabel('Drop Magnitude (Hz)')
title('τ vs B')

% Print statistical relationships
fprintf('\nRelationships with Model Parameters:\n')
[r, p] = corrcoef(log10(brsty + eps), modelParams(:,3));
fprintf('Burstiness vs Time Constant: r=%.2f, p=%.3f\n', r(1,2), p(1,2))

[r, p] = corrcoef(log10(brsty + eps), modelParams(:,2));
fprintf('Burstiness vs Drop Magnitude: r=%.2f, p=%.3f\n', r(1,2), p(1,2))

[r, p] = corrcoef(modelParams(:,3), modelParams(:,2));
fprintf('Time Constant vs Drop Magnitude: r=%.2f, p=%.3f\n', r(1,2), p(1,2))

fprintf('\nModel Fit Quality:\n')
fprintf('Mean R²: %.2f\n', mean(modelR2))
fprintf('Median R²: %.2f\n', median(modelR2))
fprintf('Number of good fits (R² > 0.7): %d/%d\n', sum(modelR2 > 0.7), nUnits)
