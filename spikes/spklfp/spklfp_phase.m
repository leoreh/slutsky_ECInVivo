function spkPhase = spklfp_phase(sig, spkTimes, fs, varargin)
% SPKLFP_PHASE Calculates high-accuracy spkPhase coupling between Spikes and LFP.
%
%   spkPhase = SPKLFP_PHASE(sig, spkTimes, fs, varargin)
%
%   SUMMARY:
%       Quantifies the spkPhase preference of spiking activity relative to an LFP oscillation.
%       Outputs include circular statistics (Mean Angle, Vector Length) and
%       bias-corrected Pairwise Phase Consistency (PPC) metrics.
%
%   INPUTS:
%       sig         - (Vec)  Filtered LFP signal (single channel).
%       spkTimes    - (Cell) {N_units x 1} Spike times for N units [s].
%       fs          - (Num)  Sampling frequency [Hz].
%       varargin    - Parameter/Value pairs:
%           'lfpTimes'  - (Mat) [N x 2] Valid analysis intervals [s].
%                         *Critical for PPC1/PPC2*: Defines "trials" for
%                         between-trial consistency (minimizing burst bias).
%           'nPerms'    - (Num) Number of permutations for statistical testing. {1000}.
%           'minSpks'   - (Num) Min spikes required to analyze a unit. {10}.
%           'sigOffset' - (Num) Time offset of signal start [s]. {0}.
%           'verbose'   - (Log) Print progress? {true}.
%
%   OUTPUTS:
%       spkPhase       - (Struct) Results with [1 x N_units] fields:
%           .theta      - Mean spkPhase angle [rad].
%           .mrl        - Mean Resultant Length (0-1).
%           .ppc0       - Standard PPC (All-to-All). Unbiased by N, biased by bursts.
%           .ppc1       - Between-Trial PPC. Unbiased by N AND bursts. (Req lfpTimes).
%           .ppc2       - Weighted Between-Trial PPC. (Req lfpTimes).
%           .kappa      - Von Mises concentration parameter.
%           .pval       - Empirical p-value from permutation test.
%           .zScore     - Z-score of MRL relative to permutation distribution.
%           .spkPh      - Cell array of phases for each valid spike.
%
%   DEPENDENCIES:
%       None.
%
%   REFERENCES:
%       Vinck et al. (2012) "Improved measures of spkPhase-coupling..."
%
%   HISTORY:
%       Updated: 23 Jan 2026

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'sig', @isnumeric);
addRequired(p, 'spkTimes', @iscell);
addRequired(p, 'fs', @isnumeric);
addParameter(p, 'lfpTimes', [], @isnumeric);
addParameter(p, 'nPerms', 0, @isnumeric);
addParameter(p, 'minSpks', 10, @isnumeric);
addParameter(p, 'sigOffset', 0, @isnumeric);

parse(p, sig, spkTimes, fs, varargin{:});

lfpTimes = p.Results.lfpTimes;
nPerms = p.Results.nPerms;
minSpks = p.Results.minSpks;
sigOffset = p.Results.sigOffset;

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
spkPhase = struct();
spkPhase.theta = nan(nUnits, 1);
spkPhase.mrl = nan(nUnits, 1);
spkPhase.ppc0 = nan(nUnits, 1);
spkPhase.ppc1 = nan(nUnits, 1);
spkPhase.ppc2 = nan(nUnits, 1);
spkPhase.kappa = nan(nUnits, 1);
spkPhase.pval = nan(nUnits, 1);
spkPhase.zScore = nan(nUnits, 1);
spkPhase.nSpks = zeros(nUnits, 1);
spkPhase.spkPh = cell(nUnits, 1);


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
    spkPhase.nSpks(iUnit) = nSpks;

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

    spkPhase.spkPh{iUnit} = angle(zSpks);

    % 5. Observed Stats (MRL & PPC)

    % A. Mean Resultant Length (MRL) and Angle
    grandSum = sum(zSpks);
    R = abs(grandSum) / nSpks;
    spkPhase.mrl(iUnit) = R;
    spkPhase.theta(iUnit) = mod(angle(grandSum), 2*pi);

    % B. PPC Calculation
    % We need trial IDs for each spike to compute per-trial sums
    spkTrialIDs = trialIdxMap(spkIdxs);

    % Group sums by trial
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
            spkPhase.ppc0(iUnit) = (abs(grandSum)^2 - nSpks) / (nSpks * (nSpks - 1));
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
                spkPhase.ppc1(iUnit) = (S_sq - Si_sq_sum) / denom;
            end

            % --- PPC2 (Bias-Corrected Between-Trial) ---

            % Calculate Mean Vector per trial
            m_i = trialZSums ./ trialNSpks;

            S_ppc2 = sum(m_i);
            SS_ppc2 = sum(abs(m_i).^2);
            J = nActiveTrials;

            if J > 1
                spkPhase.ppc2(iUnit) = (abs(S_ppc2)^2 - SS_ppc2) / (J * (J - 1));
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
    spkPhase.kappa(iUnit) = k;

    % 6. Permutation Test
    % Based on MRL (equivalent significance to PPC for testing H0)

    if nPerms > 0
        randIdx = randi(length(permPool), nSpks, nPerms);
        randPh = permPool(randIdx);

        zSumP = sum(exp(1i * randPh), 1);
        permMRL = abs(zSumP) / nSpks;

        spkPhase.pval(iUnit) = mean(permMRL >= R);

        muPerm = mean(permMRL);
        sigPerm = std(permMRL);
        spkPhase.zScore(iUnit) = (R - muPerm) / sigPerm;
    end
end

end         % EOF
