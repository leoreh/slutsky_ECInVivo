function lmeTbl = lme_mdl2tbls(lmeMdl)
% ORGANIZE_LME_FOR_EXCEL Extracts key information from an LME model
% into a single wide table suitable for Excel export.
%
% INPUT:
%   lmeMdl - A LinearMixedModel object (from fitlme).
%
% OUTPUT:
%   outputTable - A MATLAB table with combined LME model information.

% Extract Model-Level Information 
formulaStr = lmeMdl.Formula.char;
fitMet = lmeMdl.FitMethod;
numObs = lmeMdl.NumObservations;

% Model fit params
fitCrit = lmeMdl.ModelCriterion;
lmeTbl.Fit = dataset2table(fitCrit);

% Fixed effects
lmeTbl.Fixed = dataset2table(lmeMdl.Coefficients);

% Covariance Parameters and Error
[~, ~, resStats] = covarianceParameters(lmeMdl);
lmeTbl.Cov = dataset2table(resStats{1});
lmeTbl.Err = dataset2table(resStats{2});

% Random effects 
[~, ~, randStats] = randomEffects(lmeMdl, 'DFMethod', 'Residual');
lmeTbl.Rand = dataset2table(randStats);


end