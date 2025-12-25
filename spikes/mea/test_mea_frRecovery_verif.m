%% Test mea_frRecovery Verification
% Creates synthetic data and tests the new mea_frRecovery function with flgPlot

clear; clc; close all;

%% 1. Generate Synthetic Data
nUnits = 10;
nTime = 300; % 300 minutes (5 hours)
t = linspace(-60, 240, nTime); % Minutes
t = t * 60; % Convert to seconds for input
idxPert = find(t>=0, 1);

fr = zeros(nUnits, nTime);

% Generate some patterns
for i = 1:nUnits
    % Baseline
    fr(i, 1:idxPert-1) = 5 + randn(1, idxPert-1)*0.5;

    % Drop at perturbation
    fr(i, idxPert:idxPert+10) = linspace(5, 0.5, 11);

    % Recovery
    if mod(i, 2) == 0
        % Even units recover
        tRec = 1:(nTime - (idxPert+10));
        fr(i, idxPert+11:end) = 0.5 + (5-0.5)*(1-exp(-0.02 * tRec));
    else
        % Odd units don't recover fully
        fr(i, idxPert+11:end) = 0.5;
    end

    % Add noise
    fr(i, :) = fr(i, :) + randn(1, nTime) * 0.2;
    fr(i, fr(i,:)<0) = 0;
end

%% 2. Run mea_frRecovery
fprintf('Running mea_frRecovery with flgPlot=true...\n');

try
    idxTrough = idxPert + 5; % Mock trough index
    rcv = mea_frRecovery(t, fr, 'idxTrough', idxTrough, ...
        'flgPlot', true, 'binSize', 60);

    fprintf('Function executed successfully.\n');
catch ME
    fprintf('Error: %s\n', ME.message);
    rethrow(ME);
end

%% 3. Verify Output
assert(isfield(rcv, 'frBsl'), 'Missing frBsl');
assert(isfield(rcv, 'pertDepth'), 'Missing pertDepth');
assert(isfield(rcv, 'spkDfct'), 'Missing spkDfct');

fprintf('Output structure verification passed.\n');
fprintf('Please inspect the generated figure for 5 unit plots.\n');
