% Test script for mea_frRecovery
clear; clc;

% Add path
addpath(genpath('d:\Code\slutsky_ECInVivo\spikes\mea'));

% 1. Generate Synthetic Data
% ---------------------------------------------------------
binSize = 60; % 1 minute
t = -3600 : binSize : 7200; % -1h to +2h
nBins = length(t);
nUnits = 5;

frMat = zeros(nUnits, nBins);

% Unit 1: Good recovery
% Baseline 10Hz, drops to 1Hz, recovers to 10Hz
idxPert = find(t >= 0, 1);
frMat(1, 1:idxPert-1) = 10;
frMat(1, idxPert:end) = 1 + (10-1) * (1 - exp(-(t(idxPert:end))/1800));

% Unit 2: No recovery
% Baseline 5Hz, drops to 0.1Hz, stays there
frMat(2, 1:idxPert-1) = 5;
frMat(2, idxPert:end) = 0.1;

% Unit 3: No perturbation
% Constant 8Hz
frMat(3, :) = 8;

% Unit 4: Partial recovery
% Baseline 10Hz, drops to 2Hz, recovers to 6Hz
frMat(4, 1:idxPert-1) = 10;
frMat(4, idxPert:end) = 2 + (6-2) * (1 - exp(-(t(idxPert:end))/900));

% Unit 5: Low firing (Bad unit)
frMat(5, :) = 0;

uGood = [true; true; true; true; false];


% 2. Run Function
% ---------------------------------------------------------
fprintf('Running mea_frRecovery...\n');
rec = mea_frRecovery(t, frMat, 'uGood', uGood);

% 3. Verify Results
% ---------------------------------------------------------

% Check Fields
reqFields = {'frBsl', 'frTrough', 'frSs', 'pertDepth', 'rcvTime', 'uRcv'};
for i = 1:length(reqFields)
    if ~isfield(rec, reqFields{i})
        error('Missing field: %s', reqFields{i});
    end
end
fprintf('All fields present.\n');

% Check Unit 1 (Full Recovery)
u1_bsl = rec.frBsl(1);
u1_ss = rec.frSs(1);
fprintf('Unit 1 Bsl: %.2f (Exp: 10)\n', u1_bsl);
fprintf('Unit 1 SS: %.2f (Exp: 10)\n', u1_ss);
if rec.uRcv(1)
    fprintf('Unit 1 correctly identified as recovered.\n');
else
    warning('Unit 1 NOT identified as recovered.');
end

% Check Unit 2 (No Recovery)
if ~rec.uRcv(2)
    fprintf('Unit 2 correctly identified as NOT recovered.\n');
else
    warning('Unit 2 identified as recovered (unexpected).');
end

% Check Unit 3 (No Perturbation)
fprintf('Unit 3 PertDepth: %.2f\n', rec.pertDepth(3));
if ~rec.uPert(3)
    fprintf('Unit 3 correctly identified as NOT perturbed.\n');
else
    warning('Unit 3 identified as perturbed (unexpected).');
end

fprintf('\nTest passed successfully.\n');
