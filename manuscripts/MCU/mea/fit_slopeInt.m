function r = fit_slopeInt(tbl, fNamePB, yVar)
% FIT_SLOPEINT fits two models per parameter set:
%
%   1. RAW model (flgStnd = false):
%      yVar ~ pBurst * genotype + (1|sbjID)
%      Intercept = genotype effect at pBurst = 0 (fixed biological ref).
%      Slope = pBurst x genotype (unstandardized, scale-dependent).
%
%   2. STANDARDIZED model (flgStnd = true):
%      Same formula, but lme_analyse Z-scores all continuous predictors
%      via tbl_trans. Slope is in SD units (comparable across definitions).
%
% FR is excluded from both models. Rationale: FR x genotype is
% non-significant, and at wide burst definitions pBurst becomes collinear
% with FR. Excluding FR allows pBurst to absorb activity-dependent variance
% at wide thresholds — which is the point of the analysis (testing whether
% the intercept approaches zero when pBurst captures all activity).
%
% RETURNS:
%   r.intEst    / r.intP      — intercept from raw model
%   r.slopeEst  / r.slopeP    — slope from raw model
%   r.slopeEstZ / r.slopePZ   — slope from standardized model
%   NaN if the model fails or coefficients are not found.

    r = struct('intEst', NaN, 'intP', NaN, ...
               'slopeEst', NaN, 'slopeP', NaN, ...
               'slopeEstZ', NaN, 'slopePZ', NaN);

    % Check variance
    pBraw = tbl.(fNamePB);
    if std(pBraw, 'omitnan') == 0 || all(isnan(pBraw))
        return
    end

    frml = sprintf('%s ~ %s * genotype + (1|sbjID)', yVar, fNamePB);

    % --- Model 1: RAW (intercept at pBurst = 0) ---
    try
        [mdl, ~, ~] = lme_analyse(tbl, frml, ...
            'dist', 'normal', 'flgStnd', false, ...
            'flgPlot', false, 'verbose', false);
    catch
        return
    end

    coefTbl  = mdl.Coefficients;
    idxSlope = find(contains(coefTbl.Name, fNamePB) & contains(coefTbl.Name, ':'));
    idxGeno  = find(contains(coefTbl.Name, 'genotype') & ~contains(coefTbl.Name, ':'));

    if ~isempty(idxSlope) && ~isempty(idxGeno)
        r.intEst   = coefTbl.Estimate(idxGeno(1));
        r.intP     = coefTbl.pValue(idxGeno(1));
        r.slopeEst = coefTbl.Estimate(idxSlope(1));
        r.slopeP   = coefTbl.pValue(idxSlope(1));
    end

    % --- Model 2: STANDARDIZED (slope in SD units) ---
    try
        [mdlZ, ~, ~] = lme_analyse(tbl, frml, ...
            'dist', 'normal', 'flgStnd', true, ...
            'flgPlot', false, 'verbose', false);
    catch
        return
    end

    coefTblZ  = mdlZ.Coefficients;
    idxSlopeZ = find(contains(coefTblZ.Name, fNamePB) & contains(coefTblZ.Name, ':'));

    if ~isempty(idxSlopeZ)
        r.slopeEstZ = coefTblZ.Estimate(idxSlopeZ(1));
        r.slopePZ   = coefTblZ.pValue(idxSlopeZ(1));
    end
end
