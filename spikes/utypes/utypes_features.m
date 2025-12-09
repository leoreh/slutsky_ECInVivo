function fetTbl = utypes_features(varargin)

% UTYPES_FEATURES calculates waveform features to classify neuronal units.
% It computes metrics, distributions, and basic classification for Regular Spiking (RS)
% vs Fast Spiking (FS) units, creating a table of features for further analysis.
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders. If empty,
%                uses current directory.
%   flgPlot      (logical) Flag to generate visualization plots
%                {true}
%   v            (struct) Pre-loaded data structure.
%
% OUTPUT:
%   fetTbl       (table) Table containing data for each unit:
%                - UnitID: Unique unit identifier
%                - initType: Classification (1=RS, 2=FS) based on TP duration
%                - bad: Boolean flag for bad units (inverted or nan tp)
%                - Feature columns...

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addOptional(p, 'basepaths', {}, @(x) iscell(x));
addOptional(p, 'flgPlot', true, @islogical);
addOptional(p, 'v', struct(), @isstruct);

parse(p, varargin{:});
basepaths = p.Results.basepaths;
flgPlot = p.Results.flgPlot;
v = p.Results.v;

%% ========================================================================
%  CREAT TABLE
%  ========================================================================

% load state vars and concatenate metrics
if isempty(fieldnames(v))
    if isempty(basepaths)
        error('Either basepaths or v must be provided.');
    end
    vars = {'swv_metrics', 'st_metrics', 'fr', 'units'};
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);
end

% Get mouse names for table
tagFiles = struct();
if ~isempty(basepaths)
    tagFiles.Mouse = get_mname(basepaths);
    [~, fileNames] = fileparts(basepaths);
    tagFiles.File = fileNames;
elseif isfield(v, 'name')
    warning('Basepaths not provided; Mouse names might be missing in table.');
end

% Define variable map for v2tbl
% Spike times features
varMap = struct();
varMap.mfr = 'fr.mfr';
varMap.royer = 'st.royer';
varMap.royer2 = 'st.royer2';
varMap.lidor = 'st.lidor';
varMap.mizuseki = 'st.mizuseki';

% Waveform features
varMap.tp = 'swv.tp';
varMap.inverted = 'swv.inverted';
varMap.tpAmp = 'swv.tpAmp';
varMap.tpRatio = 'swv.tpRatio';
varMap.tpSlope = 'swv.tpSlope';
varMap.spkw = 'swv.spkw';
varMap.asym = 'swv.asym';
varMap.hpk = 'swv.hpk';
varMap.tailSlope = 'swv.tailSlope';
varMap.tailAmp = 'swv.tailAmp';

% Create table using v2tbl
fetTbl = v2tbl('v', v, 'varMap', varMap, 'tagFiles', tagFiles);

%% ========================================================================
%  PRELIMINARY CLASSIFICATION
%  ========================================================================

% Identify bad units (inverted or missing tp)
fetTbl.bad = fetTbl.inverted | isnan(fetTbl.tp);
idxGood = ~fetTbl.bad;

if sum(idxGood) == 0
    warning('No good units found (non-inverted).');
    return;
end

% For non-inverted spikes, perform classification using trough-to-peak
% RS: tp >= 0.4 ms (Regular Spiking, typically Pyramidal)
% FS: tp < 0.4 ms (Fast Spiking, typically Interneuron)
initType = zeros(height(fetTbl), 1);
initType(idxGood & (fetTbl.tp >= 0.4)) = 1;  % RS
initType(idxGood & (fetTbl.tp < 0.4)) = 2;   % FS

% Create a categorical variable for plotting and storage
fetTbl.initType = categorical(initType, [0 1 2], {'Other', 'RS', 'FS'});
fetTbl.unitType = fetTbl.initType; % Alias for compatibility if needed, using initType primarily

%% ========================================================================
%  TRANSFORM SKEWED FEATURES
%  ========================================================================

% Define features to plot (numeric columns excluding tags)
pltFlds = fetTbl.Properties.VariableNames;
excludedFields = {'UnitID', 'inverted', 'Name', 'Group', 'UnitType', 'initType', 'bad'};
numericIdx = varfun(@isnumeric, fetTbl, 'OutputFormat', 'uniform');
pltFlds = pltFlds(numericIdx & ~ismember(pltFlds, excludedFields));

% Apply transformations (log-transform skewed data)
% Note: Using varsInc to limit transformation to features only
fetTbl = tbl_transform(fetTbl, 'flgLog', true, 'varsInc', pltFlds);

%% ========================================================================
%  PLOT
%  ========================================================================

if flgPlot
    % Plot histograms of waveform metrics
    figure('Name', 'Waveform Metrics Distributions', 'Position', [100 100 1200 800]);

    nFields = length(pltFlds);
    nCols = 4;
    nRows = ceil(nFields / nCols);

    tiledlayout(nRows, nCols, 'TileSpacing', 'tight', 'Padding', 'tight');

    for iFld = 1:nFields
        nexttile;
        fieldName = pltFlds{iFld};
        data = fetTbl.(fieldName)(idxGood);

        % Simple classification based metrics
        rsData = data(initType(idxGood) == 1);
        fsData = data(initType(idxGood) == 2);

        % Calculate ROC/AUC for RS vs FS separation power of this feature
        % (Just as a metric of how different they are)
        [~,~,~,auc] = perfcurve(initType(idxGood), data, 1);
        if ~isempty(rsData) && ~isempty(fsData)
            pVal = ranksum(rsData, fsData);
        else
            pVal = NaN;
        end

        histogram(data, 30, 'Normalization', 'probability');
        title(sprintf('%s\nAUC=%.3f, p=%.2e', fieldName, auc, pVal), 'Interpreter', 'none');
        grid on;
    end

    % Define colors
    % Categories in initType: {'Other', 'RS', 'FS'} -> Indices: 1, 2, 3
    % mcu_clr order: Row 1 = RS, Row 2 = FS, Row 3 = Other
    % Since only good units are plotted, no need to reorder.
    clr = mcu_clr();
    clr = clr.unitType;
    % Add alpha
    clr = [clr, repmat(0.5, 3, 1)];

    % Waveform features
    % varsInc = {'tp', 'tpAmp', 'tpRatio', 'tpSlope', 'spkw', 'asym',...
    %     'hpk', 'tailSlope', 'tailAmp'};
    % plot_corrHist(fetTbl(idxGood, :), 'varsInc', varsInc, 'grpIdx', 'UnitType', ...
    %     'clrGrp', clr);
    %
    % % Spiking Features
    % varsInc = {'mfr', 'royer', 'royer2', 'lidor', 'mizuseki'};
    % plot_corrHist(fetTbl(idxGood, :), 'varsInc', varsInc, 'grpIdx', 'UnitType', ...
    %     'clrGrp', clr);
    %
    % % Selected
    % varsRow = {'tp', 'asym', 'hpk', 'tpAmp'};
    % varsCol = {'mfr', 'lidor'};
    % plot_corrHist(fetTbl(idxGood, :), 'varsRow', varsRow, 'varsCol', varsCol,...
    %     'grpIdx', 'UnitType', 'clrGrp', clr);

    % Selected
    varsRow = {'tp', 'asym', 'hpk'};
    varsCol = {'mfr', 'lidor', 'tpAmp'};
    plot_corrHist(fetTbl(idxGood, :), 'varsRow', varsRow, 'varsCol', varsCol,...
        'grpIdx', 'unitType', 'clrGrp', clr);

end

end

% EOF

%% ========================================================================
%  NOTE: FEATURE VALIDATION VIA AUC
%  ========================================================================
%  We use the Area Under the ROC Curve (AUC) to validate whether the
%  waveform-based classification predicts other physiological properties.
%
%  1. Ground Truth:
%     The classification labels (RS vs. FS) are defined
%     solely by the Trough-to-Peak (TP) waveform duration.
%
%  2. Interpretation of AUC:
%     The AUC quantifies how well *other* metrics (e.g., Firing Rate)
%     align with these waveform labels.
%     - AUC ≈ 1.0 or 0.0: Strong Correlation. The feature is an effective
%       predictor. (Note: Values < 0.5 simply indicate an inverse
%       relationship, e.g., Higher Rate = Lower Probability of RS).
%     - AUC ≈ 0.5: Random Chance. The feature provides no discriminative
%       information for this classification.
%  ========================================================================

%% ========================================================================
%  NOTE: FEATURE SELECTION & REDUNDANCY ELIMINATION
%  ========================================================================
%  Feature selection is critical for Gaussian Mixture Model (GMM) performance.
%  GMMs rely on the Mahalanobis distance, which requires inverting the
%  covariance matrix. Highly correlated features (multicollinearity) cause
%  this matrix to become singular or unstable, leading to poor clustering.
%
%  1. Redundancy Criteria:
%     We assess redundancy using Spearman's Rank Correlation. Features with
%     |rho| > 0.85 are considered redundant.
%     - Visual Confirmation: Redundant pairs manifest as tight diagonal
%       lines in pairwise scatter plots.
%     - Independence: Effective feature pairs appear as "clouds" or
%       orthogonal distributions (rho < 0.3), providing unique information
%       to the classifier.
%
%  2. Selection Strategy (The "Survivor" Rule):
%     When a redundant pair is identified, we retain the feature that is:
%     - More Robust: Less sensitive to noise (e.g., Prefer 'FullWidth' over
%       'TailSlope').
%     - More Gaussian: Exhibits a distribution closer to Normality (or
%       Log-Normality), fitting the GMM assumptions better.
%     - Biologically Interpretable: Prioritize metrics with direct
%       physiologial correlates (e.g., Repolarization Time).
%  ========================================================================