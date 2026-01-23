function phase = spklfp_phase(sig, spkTimes, fs, varargin)
% SPKLFP_PHASE Calculates high-accuracy phase coupling between Spikes and LFP.
%
% SUMMARY:
%   Quantifies the phase preference of spiking activity relative to an LFP oscillation.
%   Supports the Pairwise Phase Consistency (PPC) suite (Vinck et al., 2012),
%   providing unbiased estimates of phase locking robust to sample size and bursting.
%
% INPUTS:
%   sig         - (Vec) Filtered LFP signal (single channel).
%   spkTimes    - (Cell) Spike times for N units [s].
%   fs          - (Num) Sampling frequency [Hz].
%
% OPTIONAL KEY-VALUE PAIRS:
%   'lfpTimes'  - (Mat) [N x 2] Valid analysis intervals [s].
%                 **CRITICAL FOR PPC1/PPC2**: Each interval is treated as a
%                 separate "Trial". Spikes within the same interval are assumed
%                 to be dependent (bursts), while spikes across intervals are independent.
%   'nPerms'    - (Num) Number of permutations for statistical testing. {1000}
%   'minSpks'   - (Num) Minimum number of spikes required to analyze a unit. {10}
%   'sigOffset' - (Num) Time offset of the signal start relative to 0 [s]. {0}
%   'verbose'   - (Log) Print progress. {true}
%
% OUTPUT:
%   phase       - Structure with fields (1 x nUnits):
%     .theta      - Mean phase angle [rad].
%     .mrl        - Mean Resultant Length (0-1).
%     .ppc0       - Standard PPC (All-to-All). Unbiased by N, but biased by bursts.
%     .ppc1       - Between-Trial PPC. Unbiased by N AND bursts. (Requires lfpTimes).
%     .ppc2       - Weighted Between-Trial PPC. (Requires lfpTimes).
%     .kappa      - Von Mises concentration parameter.
%     .pval       - Empirical p-value from permutation test.
%     .zScore     - Z-score of MRL relative to permutation distribution.
%     .spkPh      - Cell array of phases for each valid spike.
%     .nSpks      - Number of valid spikes used.
%
% REFERENCES:
%   Vinck et al. (2012) "Improved measures of phase-coupling..."
%   FieldTrip implementation: ft_spiketriggeredspectrum_stat.m
%
% 22 Jan 2026 - Added FieldTrip-verified PPC0/1/2 logic.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'sig', @isnumeric);
addRequired(p, 'spkTimes', @iscell);
addRequired(p, 'fs', @isnumeric);
addParameter(p, 'lfpTimes', [], @isnumeric);
addParameter(p, 'nPerms', 1000, @isnumeric);
addParameter(p, 'minSpks', 10, @isnumeric);
addParameter(p, 'sigOffset', 0, @isnumeric);
addParameter(p, 'verbose', true, @islogical);

parse(p, sig, spkTimes, fs, varargin{:});

lfpTimes = p.Results.lfpTimes;
nPerms = p.Results.nPerms;
minSpks = p.Results.minSpks;
sigOffset = p.Results.sigOffset;
verbose = p.Results.verbose;

%% ========================================================================
%  PREP SIGNAL
%  ========================================================================

sig = sig(:);
nSamples = length(sig);
tStamps = (0:nSamples-1)' / fs;

% Analytic Signal (Hilbert)
% We keep the complex representation for accurate linear interpolation
sigH = hilbert(sig);

% Create Validity Mask & Trial Index Map
isValidSamp = true(nSamples, 1);
trialIdxMap = zeros(nSamples, 1); % Stores which 'lfpTime' row a sample belongs to

if ~isempty(lfpTimes)
    numTrials = size(lfpTimes, 1);

    % Adjust lfpTimes to be relative to the signal start
    relLfpTimes = lfpTimes - sigOffset;

    % Clip to signal duration
    relLfpTimes(relLfpTimes < 0) = 0;
    relLfpTimes(relLfpTimes > tStamps(end)) = tStamps(end);

    % Mask and Map
    isValidSamp = false(nSamples, 1);
    idxStart = floor(relLfpTimes(:,1) * fs) + 1;
    idxEnd = ceil(relLfpTimes(:,2) * fs) + 1;

    for iTr = 1:numTrials
        s = max(1, idxStart(iTr));
        e = min(nSamples, idxEnd(iTr));
        isValidSamp(s:e) = true;
        trialIdxMap(s:e) = iTr; % Map sample to Trial ID
    end
else
    numTrials = 1;
    trialIdxMap(:) = 1; % Treat whole recording as 1 trial if no windows provided
end

% Create a pool of "Valid Phases" for permutation testing
% We enable subsampling (10x) for efficiency without losing distribution shape
validPhases = angle(sigH(isValidSamp));
if length(validPhases) > 100000
    permPool = validPhases(1:10:end);
else
    permPool = validPhases;
end

%% ========================================================================
%  ANALYZE UNITS
%  ========================================================================

nUnits = length(spkTimes);

% Initialize Output Structure
phase = struct();
phase.theta = nan(1, nUnits);
phase.mrl = nan(1, nUnits);
phase.ppc0 = nan(1, nUnits);
phase.ppc1 = nan(1, nUnits);
phase.ppc2 = nan(1, nUnits);
phase.kappa = nan(1, nUnits);
phase.pval = nan(1, nUnits);
phase.zScore = nan(1, nUnits);
phase.nSpks = zeros(1, nUnits);
phase.spkPh = cell(1, nUnits);

if verbose
    fprintf('analyzing phase coupling for %d units (%d perms)...\n', nUnits, nPerms);
end

for iUnit = 1:nUnits

    % 1. Get Spikes
    uSpks = spkTimes{iUnit};
    if isempty(uSpks), continue; end

    % 2. Convert to Relative Time
    relSpks = uSpks - sigOffset;

    % 3. Filter Valid Spikes & Get Trial IDs
    spkIdxs = round(relSpks * fs) + 1;
    validMask = (spkIdxs >= 1) & (spkIdxs <= nSamples);

    % Apply Signal Bounds
    relSpks = relSpks(validMask);
    spkIdxs = spkIdxs(validMask);

    % Apply LFP Validity Mask
    isVal = isValidSamp(spkIdxs);
    relSpks = relSpks(isVal);
    spkIdxs = spkIdxs(isVal);

    nSpks = length(relSpks);
    phase.nSpks(iUnit) = nSpks;

    if nSpks < minSpks
        continue;
    end

    % 4. Complex Interpolation (High Accuracy)
    % Interpolate Real/Imag components separately using Linear interpolation
    rI = interp1(tStamps, real(sigH), relSpks, 'linear');
    iI = interp1(tStamps, imag(sigH), relSpks, 'linear');

    % Reconstruct Phase
    % These are complex unit vectors: z = exp(i*theta)
    zSpks = rI + 1i * iI;
    % Normalize to unit magnitude (interp1 might shrink amplitude slightly)
    zSpks = zSpks ./ abs(zSpks);

    phase.spkPh{iUnit} = angle(zSpks);

    % 5. Observed Stats (MRL & PPC)

    % A. Mean Resultant Length (MRL) and Angle
    grandSum = sum(zSpks);
    R = abs(grandSum) / nSpks;
    phase.mrl(iUnit) = R;
    phase.theta(iUnit) = mod(angle(grandSum), 2*pi);

    % B. PPC Calculation
    % We need trial IDs for each spike to compute per-trial sums
    spkTrialIDs = trialIdxMap(spkIdxs);

    % Group sums by trial
    % Use accumarray for speed: indices must be positive integers
    % spkTrialIDs contains 0? No, isVal check ensures it corresponds to isValidSamp (true),
    % where trialIdxMap is populated (>=1).
    if numTrials > 0 && ~isempty(spkTrialIDs)
        % Vectors per trial S_i
        trialZSums = accumarray(spkTrialIDs, zSpks, [numTrials 1]);
        trialNSpks = accumarray(spkTrialIDs, 1, [numTrials 1]);

        validTr = trialNSpks > 0;
        trialZSums = trialZSums(validTr);
        trialNSpks = trialNSpks(validTr);
        nActiveTrials = length(trialZSums);

        % --- PPC0 (All-to-All) ---
        % Formula: (|S|^2 - N) / (N*(N-1))
        if nSpks > 1
            phase.ppc0(iUnit) = (abs(grandSum)^2 - nSpks) / (nSpks * (nSpks - 1));
        end

        % --- PPC1 (Between-Trial) ---
        % Formula: (|S|^2 - sum(|Si|^2)) / (N^2 - sum(ni^2))
        % Matches FieldTrip ft_spiketriggeredspectrum_stat.m
        if nActiveTrials > 1
            S_sq = abs(grandSum)^2;
            Si_sq_sum = sum(abs(trialZSums).^2);

            N_sq = nSpks^2;
            ni_sq_sum = sum(trialNSpks.^2);

            denom = N_sq - ni_sq_sum;
            if denom > 0
                phase.ppc1(iUnit) = (S_sq - Si_sq_sum) / denom;
            end

            % --- PPC2 (Bias-Corrected Between-Trial) ---
            % FieldTrip Formula: (S*conj(S) - SS) / (dof * (dof-1))?
            % Wait, FieldTrip uses dof as "Trial Count" (dofS) for the weighted case?
            % Actually, looking at lines 398 in ft_spiketriggeredspectrum_stat:
            % out = (S*S' - SS) / (dofS * (dofS-1))
            % Where dofS is the SUM OF WEIGHTS (number of trials?).
            % This implies PPC2 in FieldTrip is a "Trial-Based Average".
            % Let's stick to the Vinck 2012 definition if ambiguous,
            % but FieldTrip's code is the golden reference you asked for.

            % Re-reading FieldTrip line 389:
            % For PPC2: dofS is incremented by 1 per trial.
            % So dofS = Number of Trials.
            % S   = sum( m_i ) where m_i = y_i / n_i  (Mean vector of trial i)
            % SS  = sum( m_i * m_i' )

            % Implementing FieldTrip PPC2 Logic specifically:
            % Calculate Mean Vector per trial
            m_i = trialZSums ./ trialNSpks;

            S_ppc2 = sum(m_i);
            SS_ppc2 = sum(abs(m_i).^2);
            J = nActiveTrials;

            if J > 1
                phase.ppc2(iUnit) = (abs(S_ppc2)^2 - SS_ppc2) / (J * (J - 1));
            end
        end
    end

    % Kappa (Von Mises Concentration)
    if R < 0.53
        k = 2 * R + R^3 + (5 * R^5) / 6;
    elseif R < 0.85
        k = -0.4 + 1.39 * R + 0.43 / (1 - R);
    else
        k = 1 / (R^3 - 4 * R^2 + 3 * R);
    end
    phase.kappa(iUnit) = k;

    % 6. Permutation Test
    % Based on MRL (equivalent significance to PPC for testing H0)

    if nPerms > 0
        randIdx = randi(length(permPool), nSpks, nPerms);
        randPh = permPool(randIdx);

        zSumP = sum(exp(1i * randPh), 1);
        permMRL = abs(zSumP) / nSpks;

        phase.pval(iUnit) = mean(permMRL >= R);

        muPerm = mean(permMRL);
        sigPerm = std(permMRL);
        phase.zScore(iUnit) = (R - muPerm) / sigPerm;
    end
end

end         % EOF
