function fetTbl = utypes_features(varargin)
% UTYPES_FEATURES Analyzes and visualizes waveform features for neuronal units.
%
% SUMMARY:
% This function calculates and visualizes various waveform features and their
% relationships for neuronal units. It computes metrics like AUC, p-values,
% and covariances between features, and generates visualization plots.
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders. If empty,
%                uses current directory.
%   flgPlot      (logical) Flag to generate visualization plots
%                {true}
%
% OUTPUT:
%   fetTbl       (table) Table containing metrics for each feature and feature pair:
%                - Name: Feature or pair name
%                - AUC: Area under ROC curve
%                - PValue: Statistical significance
%                - CovPyr: Covariance for putative pyramidal cells
%                - CovInt: Covariance for putative interneurons
%
% DEPENDENCIES:
%   basepaths2vars, catfields
%
% HISTORY:
%   Aug 2024 (AI Assisted) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments using inputParser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepaths', @(x) iscell(x));
addOptional(p, 'flgPlot', true, @islogical);

parse(p, varargin{:});
basepaths = p.Results.basepaths;
flgPlot = p.Results.flgPlot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load state vars and concatenate metrics
vars = {'swv_metrics', 'st_metrics', 'fr', 'units'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

% concatenate metrics across recordings
swv = catfields([v.swv], 2);  % waveform metrics
st = catfields([v.st], 2);    % spike train metrics
fr = catfields([v.fr], 1);    % firing rate metrics

% get number of units
nUnits = length(swv.tp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial classification using threshold (for feature analysis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Identify non-inverted spikes
idxGood = ~swv.inverted & ~isnan(swv.tp);  
nGood = sum(idxGood);

% For non-inverted spikes, perform classification using trough-to-peak
% pPYR: tp >= 0.4, pINT: tp < 0.4
unitType = zeros(1, nGood);  % Initialize all as 0 (inverted/unknown)
unitType(swv.tp(idxGood) >= 0.4) = 1;  % pPYR
unitType(swv.tp(idxGood) < 0.4) = 2;   % pINT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate metrics for all features and pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get all field names from swv struct
swvFlds = fieldnames(swv);

% Filter fields to only include row vectors
pltFlds = {};
for iFld = 1:length(swvFlds)
    fieldName = swvFlds{iFld};
    data = swv.(fieldName);
    if isrow(data) && isnumeric(data)
        pltFlds{end + 1} = fieldName;
    end
end

% Define pairs of features to plot (indices into pltFlds)
featurePairs = [
    1, 2;    
    1, 3;    
    1, 4;    
    1, 5;    
    1, 6;    
    1, 7; 
    1, 8;
    1, 9;  
    1, 12;
    6, 7;
    2, 4;
];

% Initialize arrays for feature metrics
nFields = length(pltFlds);
aucVals = zeros(nFields, 1);
pVals = zeros(nFields, 1);

% Calculate metrics for each feature
for iFld = 1:nFields
    data = swv.(pltFlds{iFld})(idxGood);
    
    % Calculate AUC
    [~, ~, ~, aucVals(iFld)] = perfcurve(unitType, data, 1);
    
    % Calculate p-value
    pyrData = data(unitType == 1);
    intData = data(unitType == 2);
    [pVals(iFld), ~] = ranksum(pyrData, intData);
end

% Create feature metrics table
featureTable = table(pltFlds', aucVals, pVals, ...
    'VariableNames', {'Feature', 'AUC', 'PValue'});

% Initialize arrays for pair metrics
nPairs = size(featurePairs, 1);
pairNames = cell(nPairs, 1);
covPyrVals = zeros(nPairs, 1);
covIntVals = zeros(nPairs, 1);

% Calculate metrics for each pair
for iPair = 1:nPairs
    fld1 = pltFlds{featurePairs(iPair, 1)};
    fld2 = pltFlds{featurePairs(iPair, 2)};
    pairNames{iPair} = sprintf('%s_vs_%s', fld1, fld2);
    
    % Get data for both features
    data1 = swv.(fld1)(idxGood);
    data2 = swv.(fld2)(idxGood);
    
    % Calculate covariance for each class
    pyrIdx = unitType == 1;
    intIdx = unitType == 2;
    
    covPyr = cov(data1(pyrIdx), data2(pyrIdx), 'omitrows');
    covInt = cov(data1(intIdx), data2(intIdx), 'omitrows');
    
    covPyrVals(iPair) = covPyr(1,2);
    covIntVals(iPair) = covInt(1,2);
end

% Create pair metrics table
pairTable = table(pairNames, covPyrVals, covIntVals, ...
    'VariableNames', {'Pair', 'CovPyr', 'CovInt'});

% Combine tables into a single metrics table
fetTbl = table('Size', [nFields + nPairs, 5], ...
    'VariableTypes', {'cell', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'Name', 'AUC', 'PValue', 'CovPyr', 'CovInt'});

% Fill in feature metrics
fetTbl.Name(1:nFields) = pltFlds';
fetTbl.AUC(1:nFields) = aucVals;
fetTbl.PValue(1:nFields) = pVals;
fetTbl.CovPyr(1:nFields) = NaN;
fetTbl.CovInt(1:nFields) = NaN;

% Fill in pair metrics
fetTbl.Name(nFields+1:end) = pairNames;
fetTbl.AUC(nFields+1:end) = NaN;
fetTbl.PValue(nFields+1:end) = NaN;
fetTbl.CovPyr(nFields+1:end) = covPyrVals;
fetTbl.CovInt(nFields+1:end) = covIntVals;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot histograms and feature pairs if requested
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flgPlot
    % Plot histograms of waveform metrics
    figure('Name', 'Waveform Metrics Distributions', 'Position', [100 100 1200 800]);
    
    % Calculate layout
    nCols = 4;
    nRows = ceil(nFields / nCols);
    
    % Create tiled layout
    hTile = tiledlayout(nRows, nCols, 'TileSpacing', 'tight', 'Padding', 'tight');
    
    % Plot histogram for each field
    for iFld = 1:nFields
        nexttile;
        
        % Get data for this field
        data = swv.(pltFlds{iFld})(idxGood);
        
        % Get stored metrics
        auc = fetTbl.AUC(iFld);
        pval = fetTbl.PValue(iFld);
        
        % Create histogram
        histogram(data, 30, 'Normalization', 'probability');
        
        % Create title with AUC and p-value
        titleStr = sprintf('%s\nAUC=%.3f, p=%.2e', pltFlds{iFld}, auc, pval);
        title(titleStr, 'Interpreter', 'none');
        xlabel('Value');
        ylabel('Probability');
        grid on;
    end
    
    % Plot feature pairs using gscatter
    figure('Name', 'Feature Pair Relationships', 'Position', [100 100 1200 800]);
    
    % Calculate layout
    nCols = 3;
    nRows = ceil(nPairs / nCols);
    
    % Create tiled layout
    hTile = tiledlayout(nRows, nCols, 'TileSpacing', 'tight', 'Padding', 'tight');
    
    % Plot each pair
    for iPair = 1:nPairs
        nexttile;
        
        % Get feature names    
        fld1 = pltFlds{featurePairs(iPair, 1)};
        fld2 = pltFlds{featurePairs(iPair, 2)};
        
        % Get data for both features
        data1 = swv.(fld1)(idxGood);
        data2 = swv.(fld2)(idxGood);
        
        % Get stored metrics
        pairIdx = nFields + iPair;
        covPyr = fetTbl.CovPyr(pairIdx);
        covInt = fetTbl.CovInt(pairIdx);
        
        % Create scatter plot colored by unitType
        gscatter(data1, data2, unitType, [0 0 1; 1 0 0], '.', [], 'off');
        
        % Add labels and title with covariance information
        xlabel(fld1, 'Interpreter', 'none');
        ylabel(fld2, 'Interpreter', 'none');
        titleStr = sprintf('%s vs %s\nCov(pPYR)=%.3f, Cov(pINT)=%.3f', ...
            fld1, fld2, covPyr, covInt);
        title(titleStr, 'Interpreter', 'none');
        grid on;
        
        % Add legend only to first plot
        if iPair == 1
            legend({'pPYR', 'pINT'}, 'Location', 'best');
        end
    end
end

end

% EOF 