function mdl = lme_fit(tbl, frml, varargin)
% LME_FIT Fits a Linear or Generalized Linear Mixed-Effects Model.
%
%   MDL = LME_FIT(TBL, FRML, ...) fits a mixed-effects model to the data in
%   table TBL using the formula FRML. Automatically selects between fitlme
%   and fitglme based on the specified distribution.
%
%   INPUTS:
%       tbl         - (table) Data table containing response and predictors.
%       frml        - (char/string) Model formula (Wilkinson notation).
%       varargin    - (param/value) Optional parameters:
%                     'dist'      : (char) Distribution of response {'Normal'},
%                                   'Binomial', 'Poisson', 'Gamma', 'InverseGaussian'.
%                     'link'      : (char) Link function. Auto-selected if empty.
%                     'fitMethod' : (char) Fitting method.
%                                   'REML' (default for LME), 'ML'.
%                                   'Laplace', REMPL (default for GLME).
%
%   OUTPUTS:
%       mdl         - (object) Fitted GeneralizedLinearMixedModel or LinearMixedModel.
%
%   See also: FITLME, FITGLME, LME_EFFECTS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'frml', @(x) ischar(x) || isstring(x));
addParameter(p, 'dist', 'Normal', @ischar);
addParameter(p, 'link', '', @ischar);
addParameter(p, 'fitMethod', '', @ischar);

parse(p, tbl, frml, varargin{:});
dist = p.Results.dist;
link = p.Results.link;
fitMethod = p.Results.fitMethod;


%% ========================================================================
%  FIT MODEL
%  ========================================================================

% Determine if GLME is needed
flgG = ~strcmpi(dist, 'Normal');

if flgG
    % --- Generalized Linear Mixed-Effects Model ---

    % Auto-select link if not provided
    if isempty(link)
        switch lower(dist)
            case 'binomial'
                link = 'Logit';
            case {'poisson', 'gamma', 'inversegaussian'}
                link = 'Log';
            otherwise
                link = 'Identity';
        end
    end

    % Default FitMethod for GLME
    if isempty(fitMethod)
        fitMethod = 'REMPL';
    end

    % Attempt 1: Standard Fit
    try
        % Temporarily treat convergence warning as error to catch
        warnId = 'stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:Message_PLUnableToConverge';
        warnState = warning('error', warnId);
        cleanupObj = onCleanup(@() warning(warnState));

        mdl = fitglme(tbl, frml, ...
            'Distribution', dist, ...
            'Link', link, ...
            'FitMethod', fitMethod);

    catch ME
        % If convergence failure and using PL method, retry with Laplace
        if contains(ME.identifier, warnId) && any(strcmpi(fitMethod, {'REMPL', 'MPL'}))
            
            warning('LME_FIT: PL failed to converge. Retrying with Laplace.');
            
            % Attempt 2: Relaxed Fit (Switching to Laplace)
            mdl = fitglme(tbl, frml, ...
                'Distribution', dist, ...
                'Link', link, ...
                'FitMethod', 'Laplace', ...
                'PLIterations', 200, ...
                'PLTolerance', 1e-6);
        else
            rethrow(ME);
        end
    end


else
    % --- Linear Mixed-Effects Model ---

    % Default FitMethod for LME
    if isempty(fitMethod)
        fitMethod = 'REML';
    end

    mdl = fitlme(tbl, frml, ...
        'FitMethod', fitMethod);
end

end         % EOF




