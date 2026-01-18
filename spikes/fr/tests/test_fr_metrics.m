function test_fr_metrics()
% EXPERIMENT: Verify NaN handling and Metric Correctness
%
%   1. Creates synthetic spike data for:
%      - Active correlated neurons
%      - Anti-correlated neurons
%      - Silent neuron
%      - Constant firing neuron (The bug source)
%
%   2. Runs fr_corr and fr_network
%
%   3. Checks:
%      - Silent/Constant units are excluded (NaN in output)
%      - mccAbs > mcc (due to anti-correlation)
%      - funcon is populated in frNet

rng(42); % Reproducibility

%% 1. Generate Synthetic Data
nTime = 1000;
t = 1:nTime;

% Neuron 1 & 2: Highly Correlated + Sine wave (+ Offset)
base = sin(t/50) + 10;
y1 = base + randn(1, nTime)*0.1;
y2 = base + randn(1, nTime)*0.1;

% Neuron 3: Anti-correlated (inverted sine around +10)
y3 = (sin(t/50) * -1) + 10 + randn(1, nTime)*0.1;

% Neuron 4: Uncorrelated noise (+ Offset)
y4 = randn(1, nTime) + 10;

% Neuron 5: Constant Firing (The Bug!)
y5 = ones(1, nTime) * 10;

% Neuron 6: Silent
y6 = zeros(1, nTime);

Y = [y1; y2; y3; y4; y5; y6];

% Display Input Stats
fprintf('Input Data: 6 Neurons x %d Timepoints\n', nTime);
fprintf('Prop: [Corr, Corr, Anti, Noise, Const, Silent]\n');

%% 2. Run fr_corr
fprintf('\n--- Testing fr_corr ---\n');
cc = fr_corr(Y, 'nShuffles', 20, 'flgPlot', false);

% CHECKS
% A. Valid Units
% We expect indices 1,2,3,4 to be valid. 5 and 6 should be excluded.
expMask = logical([1 1 1 1 0 0]);
nValid = sum(~isnan(cc.funcon));

if nValid == 4
    fprintf('[PASS] Correct number of valid units (4/6). Constant/Silent removed.\n');
else
    fprintf('[FAIL] Expected 4 valid units, got %d.\n', nValid);
end

% B. Metrics
fprintf('MCC (Signed): %.4f\n', cc.mcc);
fprintf('MCC (Abs)   : %.4f\n', cc.mccAbs);

if cc.mccAbs > cc.mcc
    fprintf('[PASS] MCC Abs (%.3f) > MCC Signed (%.3f) as expected with mixed signs.\n', cc.mccAbs, cc.mcc);
else
    fprintf('[FAIL] MCC Abs should be > MCC Signed.\n');
end

%% 3. Run fr_network (Wrapper Test)
fprintf('\n--- Testing fr_network ---\n');

% Convert Y to spike times for the wrapper
spktimes = cell(6, 1);
for i = 1:6
    % simple threshold to get spikes
    spktimes{i} = find(Y(i,:) > 0.5) ./ 1000; % fake timestamps
end

% Run wrapper
frNet = fr_network(spktimes, 'binSize', 0.001, 'winLim', [0 1], 'nShuffles', 5);

% CHECK
if isfield(frNet, 'funcon') && ~all(isnan(frNet.funcon(:)))
    fprintf('[PASS] frNet.funcon exists and is populated.\n');
    disp('FunCon Output (First Chunk):');
    disp(frNet.funcon(1,:));
else
    fprintf('[FAIL] frNet.funcon is missing or all NaN.\n');
end

end
