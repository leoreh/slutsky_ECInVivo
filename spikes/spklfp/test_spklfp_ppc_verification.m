%% Test spklfp_phase.m Logic (Refined)
% This script synthetically generates LFP and Spikes to verify that PPC0,
% PPC1, and PPC2 behave as predicted by Vinck et al. (2012).
%
% REFINED: Using Low Trial Count and Tight Bursts to maximize bias visibility.

% clc; clear; close all;

%% 1. Setup Parameters
fs = 1250;
nTrials = 40; % Reduced from 200 to highlight bias
trialDur = 0.100;
fOsc = 150;
nSpksPerBurst = 8; % Increased burstiness
burstJitter = 0.0005; % Tight 0.5ms jitter (Very consistent within burst)

% Generate "Continuous" Timeline
totalTime = nTrials * (trialDur + 0.05);
t = (0:1/fs:totalTime)';
sig = sin(2*pi*fOsc*t) + 0.2*randn(size(t)); % Cleaner LFP

% Define LFP Times
lfpTimes = zeros(nTrials, 2);
cursor = 0.02;
for i = 1:nTrials
    lfpTimes(i,:) = [cursor, cursor + trialDur];
    cursor = cursor + trialDur + 0.05;
end

%% 2. Generate Spike Trains
allSpks_Null = [];
allSpks_Burst = [];
allSpks_Locked = [];

for i = 1:nTrials
    % A. NULL (Poisson)
    nRand = poissrnd(5);
    tOffset = lfpTimes(i,1) + rand(nRand,1)*trialDur;
    allSpks_Null = [allSpks_Null; tOffset];

    % B. BURST ARTIFACT (Random Phase)
    % The burst happens at a RANDOM time (Random Phase)
    % But the spikes INSIDE the burst are TIGHT (Consistent Phase relative to each other)
    tCenter = lfpTimes(i,1) + rand*trialDur;
    burst = tCenter + randn(nSpksPerBurst, 1) * burstJitter;
    allSpks_Burst = [allSpks_Burst; burst];

    % C. LOCKED (Phase Preference)
    % Pick peaks
    winIdx = t >= lfpTimes(i,1) & t <= lfpTimes(i,2);
    localT = t(winIdx);
    localSig = sig(winIdx);
    [~, locs] = findpeaks(localSig, localT);
    if ~isempty(locs)
        pks = locs(randi(length(locs), nSpksPerBurst, 1));
        sT = pks + randn(nSpksPerBurst, 1) * 0.001;
        allSpks_Locked = [allSpks_Locked; sT];
    end
end

spikes_Null{1} = allSpks_Null;
spikes_BurstArtifact{1} = allSpks_Burst;
spikes_Locked{1} = allSpks_Locked;

%% 3. Run Analysis
fprintf('Running Analysis (N=%d Trials, %d Spks/Burst)...\n', nTrials, nSpksPerBurst);

phNull = spklfp_phase(sig, spikes_Null, fs, 'lfpTimes', lfpTimes, 'verbose', false, 'nPerms', 0);
phBurst = spklfp_phase(sig, spikes_BurstArtifact, fs, 'lfpTimes', lfpTimes, 'verbose', false, 'nPerms', 0);
phLock = spklfp_phase(sig, spikes_Locked, fs, 'lfpTimes', lfpTimes, 'verbose', false, 'nPerms', 0);

%% 4. Display Results
fprintf('\n=== RESULTS ===\n');
fprintf('Metric\t\tNULL\t\tBURST(RanPhase)\tLOCKED\n');
fprintf('---------------------------------------------------\n');
fprintf('PPC0 (All)\t%.4f\t\t%.4f\t\t%.4f\n', phNull.ppc0, phBurst.ppc0, phLock.ppc0);
fprintf('PPC1 (Btwn)\t%.4f\t\t%.4f\t\t%.4f\n', phNull.ppc1, phBurst.ppc1, phLock.ppc1);
fprintf('PPC2 (Wght)\t%.4f\t\t%.4f\t\t%.4f\n', phNull.ppc2, phBurst.ppc2, phLock.ppc2);
fprintf('MRL       \t%.4f\t\t%.4f\t\t%.4f\n', phNull.mrl, phBurst.mrl, phLock.mrl);

fprintf('\nEXPECTATION:\n');
fprintf('In BURST column, PPC0 should be significantly positive (e.g. > 0.01) due to burst bias.\n');
fprintf('PPC1 should be near zero (correcting the bias).\n');
if phBurst.ppc0 > phBurst.ppc1 + 0.005
    fprintf('\nSUCCESS: PPC0 detected the burst artifact, PPC1 ignored it.\n');
else
    fprintf('\nWARNING: Difference is small. Trigger tighter bursts or fewer trials.\n');
end
