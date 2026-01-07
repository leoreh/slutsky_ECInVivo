% Verification Script for brst_detect Optimization

clear; clc;

% Add path if needed
addpath('d:\Code\slutsky_ECInVivo\spikes\brst\');

% Generate Synthetic Data
nUnits = 20;
spktimes = cell(nUnits, 1);
fprintf('Generating synthetic data...\n');
for i = 1:nUnits
    % Poisson process with some bursty injection
    rate = 5; % Hz
    T = 600; % 10 minutes
    nSpk =poissrnd(rate * T);
    st = sort(rand(nSpk, 1) * T);

    % Inject fake bursts
    nBursts = 50;
    for b = 1:nBursts
        startT = rand * (T - 1);
        % Create a burst of 5-10 spikes with small ISI
        nB = randi([5, 10]);
        bSpikes = startT + (0:nB-1)' * 0.005; % 5ms ISI
        st = [st; bSpikes];
    end
    spktimes{i} = sort(st);
end

% Parameters
params = {'isiStart', 0.015, 'isiEnd', 0.025, 'minIbi', 0.1, 'minDur', 0.02, 'minSpks', 4, 'flgPlot', false, 'flgSave', false};

% Run OLD (we need to temporarily rename or call if strictly named)
% Assuming brst_detect_old.m is available
if ~exist('brst_detect_old', 'file')
    error('brst_detect_old.m not found!');
end

fprintf('Running OLD implementation...\n');
tic;
brst_old = brst_detect_old(spktimes, params{:});
t_old = toc;
fprintf('OLD time: %.4f s\n', t_old);

fprintf('Running NEW implementation...\n');
tic;
brst_new = brst_detect(spktimes, params{:});
t_new = toc;
fprintf('NEW time: %.4f s\n', t_new);

fprintf('Speedup: %.2fx\n', t_old / t_new);

% Compare Structs
% Fields: times, nBspk, dur, freq, ibi, spktimes
% Note: spktimes might have slight ordering differences if not sorted, but they should be sorted.
% Floating point differences might exist.

fields = {'times', 'nBspk', 'dur', 'freq', 'ibi', 'spktimes'};
allMatch = true;

for i = 1:nUnits
    for f = 1:length(fields)
        fld = fields{f};
        valOld = brst_old.(fld){i};
        valNew = brst_new.(fld){i};

        if isempty(valOld) && isempty(valNew)
            continue;
        end

        % Check for numerical equality with tolerance
        if iscell(valOld)
            % Recursive check not needed for depth 1
        elseif isnumeric(valOld)
            diffVal = abs(valOld - valNew);
            if max(diffVal(:)) > 1e-9 % Tolerance for float jitter
                fprintf('Mismatch in Unit %d, Field %s\n', i, fld);
                allMatch = false;

                % Debug print
                if length(valOld) < 10
                    disp('Old:'); disp(valOld);
                    disp('New:'); disp(valNew);
                end
            end
        end
    end
end

if allMatch
    fprintf('\nSUCCESS: Outputs match exactly!\n');
else
    fprintf('\nFAILURE: Outputs do not match.\n');
end
