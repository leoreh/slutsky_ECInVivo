function unitType = utypes_classify(varargin)
% UTYPES_CLASSIFY Classifies units into putative pyramidal cells and interneurons.
%
% SUMMARY:
% This function classifies neuronal units into putative pyramidal cells (pPYR)
% and putative interneurons (pINT) using a two-step approach:
% 1. Initial classification using trough-to-peak (tp) threshold
% 2. Optional refinement using GMM clustering starting from the threshold-based
%    classification
%
% METHODOLOGY:
% Method 1 (altClassify = 1): Uses simple threshold-based classification on tp
% Method 2 (altClassify = 2): Uses GMM clustering initialized with threshold-based
%                            classification results
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders. If empty,
%                uses current directory.
%   altClassify  (numeric) Classification method to use:
%                1 - Threshold-based classification using tp
%                2 - GMM-based classification using multiple metrics, initialized
%                    with threshold-based results
%                {1}
%   flgSave      (logical) Flag to save classification results
%                {false}
%
% OUTPUT:
%   unitType     (vector) Vector of unit classifications where:
%                1 = putative pyramidal cell (pPYR)
%                2 = putative interneuron (pINT)
%
% DEPENDENCIES:
%   basepaths2vars, catfields
%
% HISTORY:
%   Aug 2024 (AI Assisted) 
%   Based in part on Oghazian et al., Front. Biomedical Tech., 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments using inputParser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepaths', @(x) iscell(x));
addOptional(p, 'altClassify', 1, @(x) ismember(x, [1, 2, 3]));
addOptional(p, 'flgSave', false, @islogical);

parse(p, varargin{:});
basepaths = p.Results.basepaths;
altClassify = p.Results.altClassify;
flgSave = p.Results.flgSave;

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
% initial classification using threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Identify non-inverted spikes
idxGood = ~swv.inverted;  
nGood = sum(idxGood);

% For non-inverted spikes, perform classification using trough-to-peak
% pPYR: tp >= 0.4, pINT: tp < 0.4
unitType = zeros(nGood, 1);  
pyrIdx = swv.tp(idxGood) >= 0.4;
unitType(pyrIdx) = 1;       % pPYR
unitType(~pyrIdx) = 2;      % pINT

% % plot bad waveforms
% wv = catfields([v.swv], 1); 
% 
% badWv = find(~idxGood);
% badWv = find(isnan(swv.rtau) & idxGood);
% hfig = figure; hAx = gca; hold on;
% for iBad = 1 : length(badWv)
%     plot(wv.wv(badWv(iBad), :))
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% refine classification if requested
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Refine classification using GMM with threshold results as starting point
switch altClassify
    case 2
        % Prepare feature matrix for GMM 
        mfr = fr.mfr(idxGood);
        fet = [swv.tp(idxGood)', st.lidor(idxGood)', mfr];
        
        % Get refined classification using GMM
        unitType = gm2units(fet, unitType);
        
    case 3
        % Prepare feature matrix for GMM 
        tailSlope = (swv.tailSlope(idxGood)');
        tailSlope = asinh(tailSlope); 
        fet = [swv.tpRatio(idxGood)', swv.asym(idxGood)',...
            swv.hpk(idxGood)', swv.tpSlope(idxGood)', swv.tp(idxGood)'];
        fet = [swv.asym(idxGood)',...
            swv.hpk(idxGood)', swv.tp(idxGood)'];
        % % fet = zscore(fet);

        % Get refined classification using GMM
        unitType = gm2units(fet, unitType);
end

% Create clean matrix for units struct
clean = false(2, nUnits);
clean(1, idxGood) = unitType == 1;  % pPYR
clean(2, idxGood) = unitType == 2;  % pINT

% Update unitType to span all units
unitType = zeros(nUnits, 1);
unitType(clean(1, :)) = 1;
unitType(clean(2, :)) = 2;

% % Plot classificaiton
% fh = figure;
% swvFld = 'tp';
% stFld = 'lidor';
% hAx = subplot(1, 2, 1); hold on
% plot_utypes('basepaths', basepaths, 'flgRaw', false,...
%     'plotType', 'scatter3', 'swvFld', swvFld, 'stFld', stFld,...
%     'unitIdx', unitType, 'hAx', hAx)
% 
% hAx = subplot(1, 2, 2); hold on
% plot_utypes('basepaths', basepaths, 'flgRaw', false,...
%     'plotType', 'wv', 'swvFld', swvFld, 'stFld', stFld,...
%     'unitIdx', unitType, 'hAx', hAx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save updated unit struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flgSave
    startIdx = 1;
    for iPath = 1 : length(basepaths)
        basepath = basepaths{iPath};
        [~, basename] = fileparts(basepath);
        uFile = fullfile(basepath, [basename, '.units.mat']);
               
        % Get number of units for this recording
        units = v(iPath).units;
        nUnits = size(units.clean, 2);
        uIdx = startIdx : startIdx + nUnits - 1;
        startIdx = uIdx(end) + 1;        

        % Add new classification field
        fieldName = sprintf('clean%d', altClassify);
        units.(fieldName) = false(2, nUnits);
        units.(fieldName)(:, :) = clean(:, uIdx);
        
        % Place altClassify3 as default
        if altClassify == 3
            units.clean = clean(:, uIdx);
        end

        % Save updated units structure
        save(uFile, 'units');
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function typeOut = gm2units(fet, typeIn)
% GM2UNITS Refines unit classification using Gaussian Mixture Model
%
% This helper function performs GMM clustering to refine the initial unit
% classification. It uses the initial classification as a starting point
% for the GMM to improve convergence and interpretability.
%
% INPUT:
%   fet      (matrix) Feature matrix for GMM [nUnits x nFeatures]
%   typeIn   (vector) Initial classification [nUnits x 1]
%                     1 = putative pyramidal cell (pPYR)
%                     2 = putative interneuron (pINT)
%
% OUTPUT:
%   typeOut  (vector) Refined classification [nUnits x 1]
%                     1 = putative pyramidal cell (pPYR)
%                     2 = putative interneuron (pINT)

% Fit GMM with 2 components using classIn as starting point
gm = fitgmdist(fet, 2, 'RegularizationValue', 0.015, 'Start', typeIn,...
    'CovarianceType', 'diagonal');

% Get cluster assignments using the cluster method
clusterLabels = cluster(gm, fet);

% Assign clusters to pPYR and pINT based on the initial classification
% We use the initial classification to determine which cluster corresponds to pPYR
% by finding which cluster has more units initially classified as pPYR
pyrCounts = zeros(1, 2);
for iUnit = 1:2
    pyrCounts(iUnit) = sum(clusterLabels == iUnit & typeIn == 1);
end
[~, pyrCluster] = max(pyrCounts);

% Create output classification vector
typeOut = ones(size(typeIn));  % Initialize all as pPYR
typeOut(clusterLabels ~= pyrCluster) = 2;  % Set non-pPYR cluster to pINT

end

% EOF 