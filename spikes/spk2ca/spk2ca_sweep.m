function results = spk2ca_sweep()
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
basepaths = basepaths(2);
vars = {'mea', 'fr'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars, 'flgPrnt', false);
nFiles = length(basepaths);

% TABLE
varMap = struct();
varMap.uGood = 'fr.uGood';
varMap.frt = 'fr.fr';
varMap.spktimes = 'mea.spktimes';
tagFiles.Name = get_mname(basepaths, 0);
tagFiles.Group = repmat(cfg.lbl.grp(1), 1, nFiles);
tagFiles.Group(contains(tagFiles.Name, 'KO')) = cfg.lbl.grp(2);
tbl = v2tbl('v', v, 'varMap', varMap, 'tagFiles', tagFiles, ...
    'uOffset', 0);
tbl(~tbl.uGood, :) = [];

% Limit spktimes to experiment
winExp = [0, 9 * 60]  * 60;
tbl.spktimes = cellfun(@(x) x(x >= winExp(1) & x <= winExp(2)), ...
    tbl.spktimes, 'UniformOutput', false);

% Time vector (1 ms resolution)
dt = 0.001;
maxTime = max(cellfun(@(x) max([0; x(:)]), tbl.spktimes));
tVec = 0:dt:maxTime;

% Time vector  (1 s resolution)
binSize  = 1;
chunks = n2chunks('n', maxTime, 'chunksize', binSize, 'lastChunk', 'exclude');
edges = [chunks(:,1)-1, chunks(:,2)];
edges = unique(edges(:))';
tBins = edges(1:end-1) + diff(edges)/2;
nBins = length(tBins);
binIdx = discretize(tVec, edges);
valid = ~isnan(binIdx);

% Time vector (1 min resolution)
% Here, FR traces are not aligned perfectly to perturbation.
xVec = v(1).fr.t / 3600;


%% ========================================================================
%  PARAMETER SWEEP
%  ========================================================================

% Parameter Spaces
Kd = [3, 6, 10];
n  = [2, 4, 8];

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
winCalc = [winBsl; length(xVec) - winBsl(2), length(xVec)];
winCalc = winCalc * 60;

% Pre-allocate
% resMat is (nUnits x 2 x nKd x nN) for parfor slicing compliance
% Dimension 2: 1=Baseline, 2=SteadyState
resMat = zeros(nUnits, 2, nKd, nN);

% Initialize parpool
if isempty(gcp('nocreate'))
    parpool('local');
end
dq = parallel.pool.DataQueue;
afterEach(dq, @(msg) fprintf('%s', msg));
fprintf('Starting Sweep [%d Units x %d Params]...\n', nUnits, nKd*nN);

parfor iUnit = 1:nUnits

    % Compute Cytosolic Calcium ONCE per Unit. Only tauC matters here.
    st = tbl.spktimes{iUnit};
    C = spk2cyto(st, 't', tVec, 'dt', dt, 'tauC', tauC, 'isiGain', isiGain);

    for iKd = 1 : nKd
        for iN = 1 : nN

            % Mito Integration
            M = cyto2mito(C, 'n', n(iN), 'Kd', Kd(iKd), ...
                'dt', dt, 'tauM', tauM, 'flgVmax', true);

            % Downsample to 1 s
            mito = accumarray(binIdx(valid)', M(valid)', [nBins 1], @mean)';

            % Mean in windows
            resMat(iUnit, 1, iKd, iN) = mean(mito(winCalc(1, 1) : winCalc(1, 2)), 'omitnan');
            resMat(iUnit, 2, iKd, iN) = mean(mito(winCalc(2, 1) : winCalc(2, 2)), 'omitnan');
        end
    end
    send(dq, sprintf('Unit %d / %d processed.\n', iUnit, nUnits));
end


% Calculate Log Ratio (Log Fold Change)

% Add floor to prevent log(0)
minVal = min(resMat(resMat > 0));
offset = minVal / 2;

% Extract matrices
mBsl = squeeze(resMat(:, 1, :, :)); % nUnits x nKd x nN
mSs  = squeeze(resMat(:, 2, :, :)); % nUnits x nKd x nN

% Log Ratio
logRatio = abs(log((mSs + offset) ./ (mBsl + offset))); % nUnits x nKd x nN

% Mean Error across Population
errMat = squeeze(mean(logRatio, 1, 'omitnan')); % Kd x N





%% ========================================================================
%  AGGREGATE & VISUALIZE
%  ========================================================================



% Find Best
[minErr, idxMin] = min(errMat(:));
[bestK_idx, bestN_idx] = ind2sub(size(errMat), idxMin);

bestKd = Kd(bestK_idx);
bestN  = n(bestN_idx);

fprintf('Optimal Parameters:\n');
fprintf('Kd = %.4f\n', bestKd);
fprintf('n  = %d\n', bestN);
fprintf('Mean Abs LogRatio = %.4f\n', minErr);

results.resMat = resMat;
results.errMat = errMat;
results.Kd     = Kd;
results.n      = n;
results.opt.Kd = bestKd;
results.opt.n  = bestN;

% Visualization
figure('Name', 'Parameter Sweep: Recovery Accuracy', 'Color', 'w');
imagesc(n, 1:nKd, errMat);
set(gca, 'YTick', 1:nKd, 'YTickLabel', arrayfun(@(x) sprintf('%.2f', x), Kd, 'UniformOutput', false));
xlabel('Hill Coeff (n)');
ylabel('Dissociation Const (Kd)');
title(sprintf('Mean LogRatio Error (Opt: Kd=%.2f, n=%d)', bestKd, bestN));
colorbar;
axis xy;

end
