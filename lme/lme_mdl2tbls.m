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
% formulaStr = lmeMdl.Formula.char;
% fitMet = lmeMdl.FitMethod;
% numObs = lmeMdl.NumObservations;

% Model fit params
fitCrit = lmeMdl.ModelCriterion;
lmeTbl.Fit = tbl_org(dataset2table(fitCrit));

% Fixed effects
lmeTbl.Fixed = tbl_org(dataset2table(lmeMdl.Coefficients));

% Covariance Parameters and Error
[~, ~, resStats] = covarianceParameters(lmeMdl);
lmeTbl.Cov = tbl_org(dataset2table(resStats{1}));
lmeTbl.Err = tbl_org(dataset2table(resStats{2}));

% Random effects 
[~, ~, randStats] = randomEffects(lmeMdl, 'DFMethod', 'Residual');
lmeTbl.Rand = tbl_org(dataset2table(randStats));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tbl = tbl_org(tbl)
% TBL_ORG Helper function to format table columns according to predefined rules

% Define columns to round to 2 decimal places
clmnRnd2 = {'Estimate', 'SE', 'SEPred', 'Statistic',...
    'tStat', 'FStat', 'Lower', 'Upper'};

% Define columns to round to 4 decimal places  
clmnRnd4 = {'pValue', 'pVal'};

% Round columns to 2 significant digits
for iCol = 1:length(clmnRnd2)
    col = clmnRnd2{iCol};
    if ismember(col, tbl.Properties.VariableNames)
        tbl.(col) = round(tbl.(col), 2, 'significant');
    end
end

% Round columns to 4 significant digits
for iCol = 1:length(clmnRnd4)
    col = clmnRnd4{iCol};
    if ismember(col, tbl.Properties.VariableNames)
        tbl.(col) = round(tbl.(col), 4, 'significant');
    end
end

end