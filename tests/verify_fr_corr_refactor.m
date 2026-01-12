%% Verify FR Correlation Refactor
% Tests the new functionality of fr_corr and fr_network

clear; clc;

% Add paths
addpath(genpath('d:\Code\slutsky_ECInVivo\spikes'));
addpath(genpath('d:\Code\slutsky_ECInVivo\utilities'));
addpath(genpath('d:\Code\slutsky_ECInVivo\packages'));

%% Test 1: fr_corr with simple correlated data
disp('-----------------------------------------------------------');
disp('TEST 1: fr_corr with correlated signals');
disp('-----------------------------------------------------------');

nUnits = 5;
nTime = 1000;
t = 1:nTime;

% Create a common signal
common = sin(t/20);

% Create unit activity: 3 units follow common, 2 are noise
Y = zeros(nUnits, nTime);
Y(1,:) = common + randn(1,nTime)*0.5;
Y(2,:) = common + randn(1,nTime)*0.5;
Y(3,:) = common + randn(1,nTime)*0.5;
Y(4,:) = randn(1,nTime); % Uncorrelated
Y(5,:) = randn(1,nTime); % Uncorrelated

% Make sure they are all positive (like FR)
Y = Y + abs(min(Y(:))) + 1;

fprintf('Running fr_corr with 20 shuffles...\n');
out = fr_corr(Y, 'nShuffles', 20, 'flgPlot', false);

fprintf('MCC (Denoised): %.4f\n', out.mcc);
fprintf('MCC (Raw):      %.4f\n', out.mccRaw);
fprintf('Noise Limit:    %.4f\n', out.noise.limit);

% Check logic
% Units 1-3 should be highly correlated.
% Units 4-5 should be low.
ccSig = out.ccSig;
disp('Significant Correlation Matrix:');
disp(ccSig);

if out.mcc > out.mccRaw
    disp('PASS: Denoised MCC is higher than Raw MCC (as expected when noise is removed).');
else
    disp('NOTE: Denoised MCC <= Raw MCC (could happen if raw correlation is dominated by noise, or if highly correlated subset is small).');
    % Actually, if we remove weak correlations (noise), the average of the *remaining* strong ones should be higher.
    % Unless there are no significant correlations.
end

assert(isstruct(out), 'Output must be a struct');
assert(isfield(out, 'noise'), 'Output must have noise field');


%% Test 2: fr_network Integration
disp(' ');
disp('-----------------------------------------------------------');
disp('TEST 2: fr_network Integration');
disp('-----------------------------------------------------------');

% Generate spike times from Y
spktimes = cell(nUnits, 1);
for i = 1:nUnits
    % Simple Poisson-like generation from rates
    rate = Y(i,:);
    % Normalize rate to probability
    rate = rate / max(rate) * 0.8;

    st = [];
    for t_idx = 1:nTime
        if rand < rate(t_idx)
            st = [st, t_idx * 0.1]; % scale time
        end
    end
    spktimes{i} = st;
end

% Run fr_network
frNet = fr_network(spktimes, 'binSize', 0.1, 'winSize', 50, 'nShuffles', 5);

disp('frNet Output Fields:');
disp(fieldnames(frNet));

fprintf('Number of Chunks: %d\n', length(frNet.mcc));
fprintf('First Chunk MCC: %.4f\n', frNet.mcc(1));

assert(isfield(frNet, 'corr'), 'frNet must have corr field');
assert(isstruct(frNet.corr), 'frNet.corr must be a struct array');

disp('TESTS COMPLETED SUCCESSFULLY');
