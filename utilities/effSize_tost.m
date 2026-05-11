function [pTost, flgEquiv, info] = effSize_tost(x1, x2, sesoi, varargin)
% EFFSIZE_TOST Two One-Sided Tests for equivalence of two groups.
%
%   [pTost, flgEquiv, info] = effSize_tost(x1, x2, sesoi, ...) tests
%   whether the difference between the means of two independent groups
%   lies inside an equivalence region (-sesoi, +sesoi).
%
%   The procedure runs two one-sided t-tests against the bounds and
%   returns the maximum of the two p-values. If pTost < alpha, the null
%   of "true effect at least as large as sesoi" is rejected and the data
%   are statistically equivalent within the bounds.
%
%   INPUT:
%       x1, x2      (vectors) data for the two groups (NaNs ignored).
%       sesoi       (double) smallest effect size of interest. Interpreted
%                   per the 'unit' option:
%                       'd'      sesoi is a Cohen's d (default).
%                       'raw'    sesoi is in the raw response units.
%
%   OPTIONAL KEY-VALUE PAIRS:
%       'unit'      'd' | 'raw' (default 'd').
%       'alpha'     (double) significance level (default 0.05).
%
%   OUTPUT:
%       pTost       max(p_lower, p_upper). Compare to alpha.
%       flgEquiv    logical, pTost < alpha (statistical equivalence).
%       info        struct with mean1, mean2, diff, se_diff, sesoi_raw,
%                   df, p_lower, p_upper, t_lower, t_upper, ci90 (the
%                   1 - 2*alpha confidence interval on the raw difference,
%                   the canonical reporting interval for TOST).
%
%   See also: EFFSIZE_D, EFFSIZE_PWR.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'x1', @isnumeric);
addRequired(p, 'x2', @isnumeric);
addRequired(p, 'sesoi', @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'unit', 'd', @(x) any(strcmpi(x, {'d', 'raw'})));
addParameter(p, 'alpha', 0.05, @(x) isnumeric(x) && x > 0 && x < 0.5);
parse(p, x1, x2, sesoi, varargin{:});

unit  = lower(p.Results.unit);
alpha = p.Results.alpha;

%% ========================================================================
%  COMPUTE
%  ========================================================================

x1 = x1(~isnan(x1));
x2 = x2(~isnan(x2));
n1 = numel(x1);
n2 = numel(x2);
df = n1 + n2 - 2;

m1 = mean(x1); m2 = mean(x2);
s1 = std(x1);  s2 = std(x2);
sPooled = sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / df);
sediff  = sPooled * sqrt(1 / n1 + 1 / n2);
diff    = m1 - m2;

% Convert sesoi to raw units if needed
if strcmp(unit, 'd')
    sesoiRaw = sesoi * sPooled;
else
    sesoiRaw = sesoi;
end

% Two one-sided t-tests
% Lower: H0 diff <= -sesoi vs H1 diff > -sesoi
tLow = (diff - (-sesoiRaw)) / sediff;
pLow = 1 - tcdf(tLow, df);

% Upper: H0 diff >= +sesoi vs H1 diff < +sesoi
tUp = (diff - sesoiRaw) / sediff;
pUp = tcdf(tUp, df);

pTost = max(pLow, pUp);
flgEquiv = pTost < alpha;

% (1 - 2*alpha) CI on raw difference (canonical TOST report)
tCrit = tinv(1 - alpha, df);
ci90 = [diff - tCrit * sediff, diff + tCrit * sediff];

%% ========================================================================
%  PACK INFO
%  ========================================================================

info = struct();
info.mean1     = m1;   info.mean2 = m2;
info.diff      = diff;
info.se_diff   = sediff;
info.sesoi_raw = sesoiRaw;
info.df        = df;
info.p_lower   = pLow; info.p_upper = pUp;
info.t_lower   = tLow; info.t_upper = tUp;
info.ci90      = ci90;

end
