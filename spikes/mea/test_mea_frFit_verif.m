%% Test mea_frFit Verification
% Creates synthetic data and tests the new mea_frFit function

clear; clc;

%% 1. Generate Synthetic Data
nUnits = 6;
nTime = 500;
t = linspace(-60, 240, nTime); % Minutes
idxPert = find(t>=0, 1);
tVec = t - t(idxPert); % Aligned time vector (0 at pert)

fr = zeros(nUnits, nTime);

% Unit 1: Perfect Exponential Recovery
fr(1, 1:idxPert-1) = 5; % Baseline
fr(1, idxPert:idxPert+20) = linspace(5, 0.5, 21); % Drop
tRec = 1:(nTime - (idxPert+20));
fr(1, idxPert+21:end) = 5 - (5-0.5)*exp(-0.02 * tRec); % Recovery

% Unit 2: No Recovery
fr(2, 1:idxPert-1) = 3;
fr(2, idxPert:end) = 0.1;

% Unit 3: Noisy Exponential
fr(3, :) = fr(1, :);
rng(1);
fr(3, :) = fr(3, :) + randn(1, nTime) * 0.2;
fr(3, fr(3,:)<0) = 0;

% Unit 4: Overshoot (Exp2-like)
fr(4, 1:idxPert-1) = 4;
fr(4, idxPert:idxPert+10) = linspace(4, 0.2, 11);
% Rise to peak
peakTime = 50;
idxPeak = idxPert + 10 + peakTime;
fr(4, idxPert+11:idxPeak) = linspace(0.2, 8, peakTime);
% Decay to SS
fr(4, idxPeak+1:end) = 8 + (4-8)*(1-exp(-0.05 * (1:(nTime-idxPeak))));

% Unit 6: Baclofen-like (Long Silence)
fr(6, 1:idxPert-1) = 6;
fr(6, idxPert:idxPert+60) = 0.01; % 60 bins of silence
tRec = 1:(nTime - (idxPert+60));
fr(6, idxPert+61:end) = 0.01 + (6-0.01)*(1-exp(-0.03 * tRec));
% Add noise
fr(6, :) = fr(6, :) + randn(1, nTime) * 0.1;
fr(6, fr(6,:)<0) = 0;


%% 2. Run mcu_detectTrough & mea_frFit
fprintf('Running mcu_detectTrough...\n');
[idxTroughVerified, ~, ~] = mcu_detectTrough(fr, idxPert, 'marginMin', 10);

fprintf('Running mea_frFit...\n');
try
    frFit = mea_frFit(fr, tVec, 'idxTrough', idxTroughVerified, ...
        'flgSave', false, 'FilterLen', 5, 'flgPlot', true);
    fprintf('Success!\n');
catch ME
    fprintf('Error: %s\n', ME.message);
    rethrow(ME);
end

%% 3. Check Dimensions
assert(size(frFit.frTrough, 1) == nUnits, 'Output units dim mismatch');
assert(size(frFit.frMdl, 1) == nUnits, 'Output model rows mismatch');
assert(size(frFit.frMdl, 2) == nTime, 'Output model time mismatch');

%% 4. Check Values
disp('Checking Unit 1 (Perfect Exp)...');
assert(frFit.goodFit(1) == true, 'Unit 1 should be good fit');
assert(abs(frFit.frBsl(1) - 5) < 0.5, 'Unit 1 Baseline incorrect');

disp('Checking Unit 2 (No Recovery)...');
% Might be bad fit mainly due to R2 if it's flat, or might be good but just very low R2
% actually if it's flat 0.1, R2 might be weird. Let's see.

disp('Checking Unit 3 (Noisy)...');
assert(frFit.goodFit(3) == true, 'Unit 3 should be good fit despite noise');

disp('Checking Unit 5 (NaNs)...');
assert(all(isnan(frFit.frMdl(5, :))), 'Unit 5 model should be NaN');

disp('Checking Unit 6 (Long Silence)...');
% We expect the fit to be decent if our new logic works, but with the OLD logic
% it usually fails to capture the delay effectively or has bad residuals.
% For now we just check it runs.
assert(frFit.goodFit(6) == true, 'Unit 6 should be good fit');


fprintf('All Basic Assertions Passed.\n');

%% 5. Plot Results
figure('Position', [100 100 1000 800]);
for i = 1:4
    subplot(2, 2, i);
    plot(tVec, fr(i, :), 'k.', 'MarkerSize', 5); hold on;
    plot(tVec, frFit.frMdl(i, :), 'r-', 'LineWidth', 2);
    title(sprintf('Unit %d: %s (R2=%.2f)', i, frFit.mdlName{i}, frFit.rsquare(i)));
    grid on;
end
