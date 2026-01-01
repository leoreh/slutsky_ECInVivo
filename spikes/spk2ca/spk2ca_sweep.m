function spk2ca_sweep()
% SPK2CA_SWEEP Performs parameter sweep (Kd, n) on real MEA data.
%
%   results = SPK2CA_SWEEP(BASEPATHS) loads spiking data from the given
%   basepaths and iterates through a range of Hill coefficients (n) and
%   dissociation constants (Kd). For each combination, it calculates the
%   Steady-State to Baseline ratio of mitochondrial calcium.
%
%   The goal is to find parameters that minimize the log-ratio (drift)
%   between baseline and steady-state periods, optimizing for 'recovery'
%   accuracy.
%
%   INPUTS:
%       basepaths   - (cell) List of recording directories.
%
%   OUTPUTS:
%       results     - (struct) Sweep results.
%                     .tbl      : (table) Unit-wise results.
%                     .sweep_Kd : (vec) Tested Kd values.
%                     .sweep_n  : (vec) Tested n values.
%                     .errMat   : (n x Kd) Population error matrix.
%
%   See also: SPK2CA, SPK2CA_SIM

%% ========================================================================
%  LOAD AND PREPARE DATA
%  ========================================================================

cfg = mcu_cfg;
basepaths = mcu_basepaths('mea_bac');
vars = {'mea', 'fr', 'frRcv'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars, 'flgPrnt', false);
nFiles = length(basepaths);

% TABLE
varMap = struct();
varMap.uGood = 'fr.uGood';
varMap.uRcv = 'rcv.uRcv';
varMap.uPert = 'rcv.uPert';
varMap.spktimes = 'mea.spktimes';
varMap.frt = 'fr.fr';
tagFiles.Name = get_mname(basepaths, 0);
tagFiles.Group = repmat(cfg.lbl.grp(1), 1, nFiles);
tagFiles.Group(contains(tagFiles.Name, 'KO')) = cfg.lbl.grp(2);
tbl = v2tbl('v', v, 'varMap', varMap, 'tagFiles', tagFiles, ...
    'uOffset', 0);
tbl.UnitID = categorical(tbl.UnitID);
tbl(~tbl.uGood, :) = [];
% tbl(~tbl.uPert, :) = [];

% Limit spktimes to experiment
winExp = [0, 9 * 60]  * 60;
spktimes = cellfun(@(x) x(x >= winExp(1) & x <= winExp(2)), ...
    tbl.spktimes, 'UniformOutput', false);

%% ========================================================================
%  TIME VECTORS
%  ========================================================================

% 1 ms resolution
dt = 0.001;
maxTime = max(cellfun(@(x) max([0; x(:)]), tbl.spktimes));
tVec = 0:dt:maxTime;

% 1 s resolution
binSize  = 1;
chunks = n2chunks('n', maxTime, 'chunksize', binSize, 'lastChunk', 'exclude');
edges = [chunks(:,1)-1, chunks(:,2)];
edges = unique(edges(:))';
tBins = edges(1:end-1) + diff(edges)/2;
nBins = length(tBins);
binIdx = discretize(tVec, edges);
valid = ~isnan(binIdx);
binIdx = binIdx(valid);

% 1 min resolution
% Here, FR traces are not aligned perfectly to perturbation.
xVec = v(1).fr.t / 3600;


%% ========================================================================
%  PARAMETER SWEEP
%  ========================================================================

% Parameter Spaces
Kd = [0.1, 1, 4, 10, 20];
n  = [2.7, 6];

nKd = length(Kd);
nN  = length(n);
nUnits = size(tbl, 1);

% Constants
isiGain  = 0;
tauC     = 0.1;
tauM     = 20;

% Windows (Indices)
idxZero = find(xVec >= 0, 1);
winBsl = [1, idxZero - 8];
winSs = [length(xVec) - winBsl(2), length(xVec) - 1];
winCalc = [winBsl; winSs] * 60 * 1000;

% Pre-allocate
% resMat is (nUnits x 2 x nKd x nN). Dimension 2: 1=Baseline, 2=SteadyState
resMat = zeros(nUnits, 2, nKd, nN);

% Initialize parpool
if isempty(gcp('nocreate'))
    parpool('local', 8);
end
dq = parallel.pool.DataQueue;
afterEach(dq, @(msg) fprintf('%s', msg));
fprintf('Starting Sweep [%d Units x %d Params]...\n', nUnits, nKd*nN);

parfor iUnit = 1:nUnits

    % Cytosolic Calcium
    st = spktimes{iUnit};
    C = spk2cyto(st, 't', tVec, 'dt', dt, 'tauC', tauC, 'isiGain', isiGain);

    tmpRes = zeros(2, nKd, nN);
    for iKd = 1 : nKd
        for iN = 1 : nN

            % Mito Integration
            M = cyto2mito(C, 'n', n(iN), 'Kd', Kd(iKd), ...
                'dt', dt, 'tauM', tauM);

            % Mean in windows
            tmpRes(1, iKd, iN) = mean(M(winCalc(1, 1) : winCalc(1, 2)), 'omitnan');
            tmpRes(2, iKd, iN) = mean(M(winCalc(2, 1) : winCalc(2, 2)), 'omitnan');
        end
    end
    resMat(iUnit, :, :, :) = tmpRes;
    send(dq, sprintf('Unit %d / %d processed.\n', iUnit, nUnits));
end


%% ========================================================================
%  ERROR
%  ========================================================================

% 1. FIRING RATE REFERENCE
% Calculate FR Log Ratio as the benchmark for homeostatic accuracy
frBsl = mean(tbl.frt(:, winBsl(1) : winBsl(2)), 2, 'omitnan');
frSs  = mean(tbl.frt(:, winSs(1) : winSs(2)), 2, 'omitnan');

% Dynamic offset for FR (1% of median baseline)
offFr = median(frBsl(frBsl > 0)) * 0.01;
logFr = abs(log((frSs + offFr) ./ (frBsl + offFr)));

% 2. MITOCHONDRIAL SWEEP ANALYSIS
errMat = zeros(nKd, nN);
riiMat = zeros(nKd, nN); % Recovery Improvement Index

for iKd = 1:nKd
    for iN = 1:nN
        % Extract data: (nUnits x 2)
        currDat = resMat(:, :, iKd, iN);
        mBsl = currDat(:, 1);
        mSs  = currDat(:, 2);

        % Define dynamic offset (1% of median baseline for this param set)
        % This prevents low-signal regimes from appearing artificially accurate
        offMito = median(mBsl(mBsl > 0), 'omitnan') * 0.01;
        if isnan(offMito) || offMito == 0, offMito = eps; end

        % SNR Mask: Only analyze units with signal significantly above the floor
        % This filters out the "vanishing signal" artifact at high Kd/n
        uValid = mBsl > (offMito * 10); 
        
        if sum(uValid) < (nUnits * 0.1)
            errMat(iKd, iN) = NaN; % Insufficient signal for reliable analysis
            continue;
        end

        % Log Ratio for Mito
        logMito = abs(log((mSs + offMito) ./ (mBsl + offMito)));

        % Recovery Improvement Index (RII)
        % Positive values mean Mito is more homeostatic than FR
        rii = logFr(uValid) - logMito(uValid);

        % Store Median Absolute Error for valid units
        errMat(iKd, iN) = median(logMito(uValid), 'omitnan');
        riiMat(iKd, iN) = median(rii, 'omitnan');
    end
end



%% ========================================================================
%  INTERACTIVE GUI
%  ========================================================================

% Flatten results into table columns for GUI
for iKd = 1 : nKd
    for iN = 1 : nN
        % Create standard variable names
        suffix = sprintf('_K%g_n%g', Kd(iKd), n(iN));
        suffix = strrep(suffix, '.', 'p'); % Safety for decimals

        tbl.(['mBsl' suffix]) = resMat(:, 1, iKd, iN);
        tbl.(['mSs' suffix])  = resMat(:, 2, iKd, iN);
    end
end

% Combine with tblu from mea_spk2ca
tblPlot = outerjoin(tbl, tblu, 'MergeKeys', true);

% Remove units with low mito
clipIdx = find(squeeze(any(any(any(resMat < 1e-20, 2), 3), 4)));
tblPlot(clipIdx, :) = [];

spk2ca_gui(tblPlot, Kd, n)

tblGUI_scatHist(tblPlot);


%% ========================================================================
%  HEATMAP: PARAMETER SWEEP RESULTS
%  ========================================================================

figure('Name', 'Parameter Sweep Results', 'Color', 'w', 'NumberTitle', 'off');

% Plot Heatmap (Transpose for n-rows, Kd-cols)
% Note: errMat is (Kd x n), so errMat' is (n x Kd) for (y, x) plotting.
valMat = errMat';
valMat = riiMat';
imagesc(Kd, n, valMat);
hCbar = colorbar;
hCbar.Label.String = 'Error (Log Ratio)';
colormap('parula');
axis xy;

% Aesthetics
title('Parameter Sweep Error');
xlabel('Dissociation Constant (Kd)');
ylabel('Hill Coefficient (n)');
set(gca, 'XTick', Kd, 'YTick', n);
set(gca, 'FontSize', 12, 'LineWidth', 1.2);

% Overlay Values
minVal = min(valMat(:));
maxVal = max(valMat(:));
rangeVal = maxVal - minVal;

for iN = 1:nN
    for iKd = 1:nKd
        val = valMat(iN, iKd);

        % Contrast Logic for Parula (Low=Dark, High=Light)
        if val < (minVal + rangeVal * 0.5)
            txtColor = 'w';
        else
            txtColor = 'k';
        end

        text(Kd(iKd), n(iN), sprintf('%.2f', val), ...
            'HorizontalAlignment', 'center', ...
            'Color', txtColor, 'FontWeight', 'bold');
    end
end


end % EOF
