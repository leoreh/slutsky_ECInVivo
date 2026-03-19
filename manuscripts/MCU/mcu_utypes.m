% mcu_cellCalss


%% ========================================================================
%  RE-ANALYZE
%  ========================================================================

% get all files in study
basepaths = mcu_basepaths('all');
nPaths = length(basepaths);

% vars
vars = {'spikes'};

% load state vars
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

for iPath = 1 : nPaths

    % files
    basepath = basepaths{iPath};
    [~, basename] = fileparts(basepath);
    cd(basepath)
    
    fr = calc_fr(spikes.times, 'basepath', basepath,...
        'graphics', false, 'binsize', 60, 'saveVar', true,...
        'smet', 'GK', 'winBL', [0, Inf], 'winCalc', [0, Inf], 'forceA', true);

    % waveform metrices
    % swv = spkwv_metrics('basepath', basepath, 'flgSave', true, 'flgForce', true);

    % Spike timing metrics
    % st = spktimes_metrics('spktimes', v(iPath).spikes.times, 'sunits', [],...
    %     'bins', {[0, Inf]}, 'flgForce', true, 'flgSave', true, 'flgAll', false);

    


end



%% ========================================================================
%  RE-CLASSIFY
%  ========================================================================

% get all files in study
basepaths = mcu_basepaths('all');
basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];

% Classify
fetSelect = {'Asym', 'Hpk', 'TP'};
rsPrior = 0.97;
regVal = 0.01;
tblUnit = utypes_classify('basepaths', basepaths(3), ...
    'fetSelect', fetSelect, 'regVal', regVal, ...
    'rsPrior', rsPrior, 'flgPlot', true);


%% ========================================================================
%  INSPECT
%  ========================================================================

basepaths = [mcu_basepaths('wt_bsl_ripp')];
[tAxis, tblUnit] = mcu_frTbl(basepaths(3), 'flgPlot', false);
utypes_gui('basepaths', basepaths(3))


% Inspect
basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];

% Grab FR vs Time data
[tAxis, tblUnit] = mcu_frTbl(basepaths, 'flgPlot', false);

% Plot classification
utypes_gui('basepaths', basepaths, 'tAxis', tAxis, 'tblUnit', tblUnit)


hFig = tblGUI_xy(tAxis, tblUnit);

% Grab to prism
idxUnits = tblUnit.UnitType == 'FS' & tblUnit.Group == 'Control';
frMat = tblUnit.FRt(idxUnits, :)';


%% ========================================================================
%  AUTOCORRELOGRAM (ACG) VISUALIZATION
%  ========================================================================

% Load unit table with narrow ACG traces.
%   acg_narrow (nunits x 201): auto-correlogram at 0.5 ms resolution,
%   computed over the full recording duration (bins = [0 Inf]).
%   xAcg (1 x 201): lag axis in milliseconds, centered at 0.
basepaths = [mcu_basepaths('wt_bsl'), mcu_basepaths('mcu_bsl')];

[tblAcg, ~, ~, xAcg] = mcu_tblVivo('basepaths', basepaths, 'presets', {'acg'}, 'flgClean', true);

% Interactive viewer: tiles = unit type, colors = genotype.
% Uses 'Spread' dispersion with arithmetic mean + SEM by default.
hFig = tblGUI_xy(xAcg.wide, tblAcg, ...
    'yVar',    'acg_narrow', ...
    'tileVar', 'unitType', ...
    'grpVar',  'Group', ...
    'xLbl',    'Lag [ms]');


%% ========================================================================
%  ACG OVER HOURS — CONVERGENCE TO DIRECT FIRING RATE
%  ========================================================================
%
% At ±500 ms lags the ACG baseline is elevated above the directly computed
% mean firing rate (N_spk / T_rec).  The non-stationarity hypothesis
% explains this through the rate-weighted sampling of reference spikes:
% because spikes are generated with probability proportional to the
% instantaneous rate lambda(t), high-rate epochs produce more reference
% spikes and the CCG baseline converges to
%
%   ACG(tau << T_fluct)  →  E[lambda^2] / E[lambda]  =  FR * (1 + CV^2)
%
% rather than FR.  The critical prediction is that ACG(tau) must decay
% back to FR_direct as tau grows beyond T_fluct, the characteristic
% timescale of firing-rate fluctuations.
%
% We test this by re-computing the ACG at hour-scale lags using the same
% normalisation as CCG (counts / N_spk / binSize), but with 60-second
% bins and lags up to 4 hours, using the FFT-based xcorr.  An edge
% correction accounts for the finite recording duration:
%
%   ACG(tau) = xcorr[lBins] * (nBins / (nBins - lBins)) / (N_spk * binSz)
%
% where the factor nBins / (nBins - lBins) corrects for the shrinking
% number of available spike pairs as tau approaches T_rec.
%
% If the long-lag ACG trace decays from the acg_wide edge value and
% converges to FR_direct, non-stationarity is the unambiguous explanation.

% Load spike times alongside st_metrics via the dedicated preset.
% flgClean is left false so all units (not just RS) are included here.
[~, ~, v_st] = mcu_tblVivo('basepaths', basepaths, 'presets', {'acg', 'spktimes'});

nPaths    = numel(basepaths);
binSz     = 60;          % [s]  coarse bin — same resolution as calc_fr
maxLag    = 4 * 3600;    % [s]  4 hours
nEdgeBins = 10;          % outermost wide-ACG bins for the edge reference

lagVec   = binSz : binSz : maxLag;   % [s]  lag axis
halfBins = maxLag / binSz;           % number of positive-lag bins

% Collector arrays — one entry per unit with a valid ACG
acgHr_all   = [];   % [nUnits x halfBins]  per-unit long-lag ACG trace  [Hz]
frDir_all   = [];   % [nUnits x 1]         N_spk / T_rec                [Hz]
acgEdge_all = [];   % [nUnits x 1]         mean of outermost acg_wide bins

for iPath = 1 : nPaths

    spikes  = v_st(iPath).spikes;
    st      = v_st(iPath).st;

    % Recording duration: last spike time across all units (spike times
    % are sorted, so this is the last element of the last non-empty cell)
    validMask = ~cellfun(@isempty, spikes.times);
    tRec = max(cellfun(@(x) x(end), spikes.times(validMask)));   % [s]

    nUnits   = numel(spikes.times);
    nLagBins = size(st.acg_wide, 2);
    edgeIdx  = [1 : nEdgeBins, nLagBins - nEdgeBins + 1 : nLagBins];

    for iUnit = 1 : nUnits

        % Skip units whose ACG was not computed (fewer than minSpkThr spikes)
        if all(isnan(st.acg_wide(iUnit, :)))
            continue
        end

        spkTimes = spikes.times{iUnit};
        nSpks    = numel(spkTimes);

        frDir   = nSpks / tRec;
        edgeVal = mean(st.acg_wide(iUnit, edgeIdx), 'omitnan');

        % Bin spike train at coarse resolution
        binEdges = 0 : binSz : (ceil(tRec / binSz) * binSz);
        N_bins   = histcounts(spkTimes, binEdges);   % spike count per 60-s bin
        nBins    = numel(N_bins);

        % Long-lag ACG via FFT-based xcorr (O(N log N)).
        %   xcorr[lBins] = sum_k N(k) * N(k + lBins) — identical in
        %   structure to the CCG pair count before normalisation.
        %   We limit the maxlag to the shorter of halfBins and nBins-1.
        nLagMax = min(halfBins, nBins - 1);
        [acgRaw, lagsBins] = xcorr(double(N_bins), nLagMax, 'none');
        acgRaw(lagsBins == 0) = 0;   % exclude zero lag (self-coincidence)

        % Positive lags only + edge correction + rate normalisation
        posIdx   = lagsBins > 0;
        lBinsVec = lagsBins(posIdx);
        corrFact = nBins ./ max(1, nBins - lBinsVec);   % finite-T correction
        acgUnit  = acgRaw(posIdx) .* corrFact / (nSpks * binSz);   % [Hz]

        % Pad with NaN for lags beyond this unit's recording duration
        acgUnit = [acgUnit, nan(1, halfBins - numel(acgUnit))];

        acgHr_all   = [acgHr_all;   acgUnit];
        frDir_all   = [frDir_all;   frDir  ];
        acgEdge_all = [acgEdge_all; edgeVal];
    end
end

% -------------------------------------------------------------------------
%  Summary
% -------------------------------------------------------------------------

acgHr_mean = mean(acgHr_all,  1, 'omitnan');
frMean     = mean(frDir_all,     'omitnan');
edgeMean   = mean(acgEdge_all,   'omitnan');

idx1h = lagVec == 3600;
fprintf('\n--- Long-lag ACG convergence ---\n')
fprintf('  ACG wide edge (tau = 0.5 s)  : %.2f Hz\n', edgeMean)
fprintf('  Long-lag ACG at tau = 1 h    : %.2f Hz\n', acgHr_mean(idx1h))
fprintf('  FR_direct = N / T_rec        : %.2f Hz\n', frMean)

% -------------------------------------------------------------------------
%  Plot: long-lag ACG population mean with reference lines
% -------------------------------------------------------------------------

hFig_hr = figure('Name', 'ACG_hours', 'Color', 'w', 'Position', [100 100 620 400]);
hAx = axes(hFig_hr);
hold(hAx, 'on')

plot(hAx, lagVec / 3600, acgHr_mean, '-', 'Color', [0.20 0.45 0.80], ...
    'LineWidth', 1.8, 'DisplayName', 'ACG(\tau)  pop. mean')

yline(hAx, edgeMean, 'r--', 'LineWidth', 1.2, ...
    'Label', sprintf('ACG_{wide} edge = %.2f Hz', edgeMean), ...
    'LabelHorizontalAlignment', 'left', 'DisplayName', 'ACG_{wide} edge')

yline(hAx, frMean, 'k--', 'LineWidth', 1.2, ...
    'Label', sprintf('FR_{direct} = %.2f Hz', frMean), ...
    'LabelHorizontalAlignment', 'left', 'DisplayName', 'FR_{direct}')

xlabel(hAx, 'Lag [h]')
ylabel(hAx, 'ACG  [Hz]')
title(hAx, 'Long-lag ACG: convergence toward FR_{direct}')
legend(hAx, 'Location', 'northeast')

