%% Demo Script for spk2ca_sim (v5: Final)
% This script demonstrates the two modes of spk2ca_sim.

clc; clear; close all;

%% 1. Mode A: Scalar TauC (Parameter Sensitivity)
% This will produce 2 figures:
% - Sensitivity: Selectivity vs Kd curves
% - Dynamics: Traces
fprintf('Running Mode A (Sensitivity Analysis)...\n');
spk2ca_sim(...
    'sweep_rate', [0.1, 1, 10], ...
    'sweep_Kd', 1:50, ...
    'tauC', 0.1, ...           % Scalar
    'flgPlot', true);

%% 2. Mode B: Vector TauC (Decay Scan)
% This will produce 3 figures:
% - Sensitivity: Selectivity vs Kd curves (Tiled)
% - Summary Analysis: Peak Sel & Gain vs TauC (Transposed)
% - Dynamics: Traces
fprintf('\nRunning Mode B (Decay Scan)...\n');
spk2ca_sim(...
    'sweep_rate', [0.1, 1, 5, 10, 20], ...
    'sweep_Kd', 1:50, ...
    'tauC', [0.02, 0.05, 0.1, 0.3, 0.5], ... % Vector
    'flgPlot', true);
