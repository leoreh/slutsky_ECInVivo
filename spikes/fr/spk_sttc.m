function [sttc, info] = spk_sttc(spk1, spk2, dt, tLim)
% SPK_STTC Calculates the Spike Time Tiling Coefficient (STTC).
%
%   [STTC, INFO] = SPK_STTC(SPK1, SPK2, DT, TLIM) computes the pairwise
%   correlation between two spike trains using the STTC metric described
%   by Cutts & Eglen (2014).
%
%   This implementation is a direct translation of the original C code
%   (spike_time_tiling_coefficient.c) to ensure exact algorithmic
%   correspondence.
%
%   INPUTS:
%       spk1    - (double) Nx1 or 1xN Spike times for train 1.
%       spk2    - (double) Mx1 or 1xM Spike times for train 2.
%       dt      - (double) Coincidence window size (delta t).
%       tLim    - (double) [Start, End] time limits of the recording.
%
%   OUTPUTS:
%       sttc    - (double) The Spike Time Tiling Coefficient (-1 to 1).
%       info    - (struct) Intermediate calculations (PA, PB, TA, TB).
%
%   REFERENCES:
%       Cutts, C. S., & Eglen, S. J. (2014). Detecting pairwise
%       correlations in spike trains: an objective comparison of methods
%       and application to optogenetics. Journal of Neuroscience, 34(43),
%       14288-14303.
%
%   See also: LME_ANALYSE

%% ========================================================================
%  INPUT PARSING
%  ========================================================================

% Ensure column vectors
spk1 = spk1(:);
spk2 = spk2(:);

N1 = length(spk1);
N2 = length(spk2);

% Sort spike trains (Crucial for the sliding window algorithm)
spk1 = sort(spk1);
spk2 = sort(spk2);

% Validate Time Limits
if nargin < 4 || isempty(tLim)
    error('[SPK_STTC] tLim [Start, End] must be provided.');
end

tStart = tLim(1);
tEnd   = tLim(2);
T      = tEnd - tStart;

% Edge Case: Empty Spike Trains
if N1 == 0 || N2 == 0
    sttc = NaN;
    info = struct('PA', NaN, 'PB', NaN, 'TA', NaN, 'TB', NaN);
    return;
end

%% ========================================================================
%  CALCULATION
%  ========================================================================

% 1. Calculate T_A and T_B
% Fraction of total time 'tiled' by spikes in each train (+/- dt)
TA_raw = run_T(N1, dt, tStart, tEnd, spk1);
TB_raw = run_T(N2, dt, tStart, tEnd, spk2);

TA = TA_raw / T;
TB = TB_raw / T;

% 2. Calculate P_A and P_B
% Proportion of spikes in one train falling within +/- dt of spikes in the other
PA_raw = run_P(N1, N2, dt, spk1, spk2);
PB_raw = run_P(N2, N1, dt, spk2, spk1);

PA = PA_raw / N1;
PB = PB_raw / N2;

% 3. Compute STTC
% Formula: 0.5 * (term1 + term2)
term1 = (PA - TB) / (1 - PA * TB);
term2 = (PB - TA) / (1 - PB * TA);

sttc = 0.5 * (term1 + term2);

% Struct for debugging/diagnostics
info = struct();
info.PA = PA;
info.PB = PB;
info.TA = TA;
info.TB = TB;

end

%% ========================================================================
%  HELPER FUNCTIONS (Match C Implementation)
%  ========================================================================

function time_A = run_T(N1, dt, start, endVal, spike_times_1)
% RUN_T Calculates the total time covered by tiles (width 2*dt)
% centered on each spike, accounting for overlaps and boundaries.

% Maximum possible time (no overlaps)
time_A = 2 * N1 * dt;

if N1 == 1
    % Singe spike case: check boundaries
    if (spike_times_1(1) - start) < dt
        time_A = time_A - start + spike_times_1(1) - dt;
    elseif (spike_times_1(1) + dt) > endVal
        time_A = time_A - spike_times_1(1) - dt + endVal;
    end
else
    % Multiple spikes: check adjacent overlaps
    % C loop: for(i=0; i<(N1-1); i++) -> 1 to N1-1 in MATLAB
    for i = 1:(N1 - 1)
        diff = spike_times_1(i+1) - spike_times_1(i);

        if diff < 2 * dt
            % Subtract overlap
            time_A = time_A - 2 * dt + diff;
        end
    end

    % Check boundaries for first and last spike
    % Note: C code handles these SEPARATELY after the loop.

    % First spike boundary
    if (spike_times_1(1) - start) < dt
        time_A = time_A - start + spike_times_1(1) - dt;
    end

    % Last spike boundary
    if (endVal - spike_times_1(N1)) < dt
        time_A = time_A - spike_times_1(N1) - dt + endVal;
    end
end
end


function Nab = run_P(N1, N2, dt, spike_times_1, spike_times_2)
% RUN_P Calculates number of spikes in train 1 that have at least
% one spike in train 2 within +/- dt.

Nab = 0;
j = 1; % Initialize j (C: j=0 -> MATLAB: 1)

for i = 1:N1
    while j <= N2
        % check every spike in train 1 to see if there's a spike in train 2
        % within dt (don't count spike pairs)

        % diff = |t1 - t2|
        if abs(spike_times_1(i) - spike_times_2(j)) <= dt
            Nab = Nab + 1;
            break; % Found a match for this i, move to next i
        elseif spike_times_2(j) > spike_times_1(i)
            % Train 2 spike is already ahead of Train 1 spike
            % Since sorted, no future T2 spike will match this T1 spike.
            break;
        else
            % spike_times_2(j) < spike_times_1[i] - dt
            % Train 2 spike is too far behind. Check next T2 spike.
            j = j + 1;
        end
    end
end
end
