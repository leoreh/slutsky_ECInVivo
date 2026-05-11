function [d, g, ci, info] = effSize_d(x1, x2, varargin)
% EFFSIZE_D Cohen's d and Hedges' g for two independent samples.
%
%   [d, g, ci, info] = effSize_d(x1, x2, ...) computes the standardized
%   mean difference between two independent groups using the pooled
%   within-group standard deviation.
%
%   INPUT:
%       x1          (vector) data for group 1 (numeric, NaNs ignored).
%       x2          (vector) data for group 2 (numeric, NaNs ignored).
%
%   OPTIONAL KEY-VALUE PAIRS:
%       'confLvl'   (double) confidence level for ci (default 0.95).
%       'flgAbs'    (logical) return |d| and |g| (default false).
%
%   OUTPUT:
%       d           Cohen's d  = (mean(x1) - mean(x2)) / sPooled.
%       g           Hedges' g  = d * (1 - 3/(4*(n1+n2)-9))  (small-n bias
%                   corrected). Recommended below n ~ 20 per group.
%       ci          [low, high] confidence interval on d at confLvl.
%                   Computed via the normal approximation
%                   se_d = sqrt((n1+n2)/(n1*n2) + d^2/(2*(n1+n2))).
%       info        struct with .mean1, .mean2, .sd1, .sd2, .sPooled,
%                   .n1, .n2, .se_d.
%
%   The convention mean(x1) - mean(x2) means a positive d indicates x1
%   exceeds x2. Effects below 0.2 are conventionally small, 0.5 medium,
%   0.8 large.
%
%   See also: EFFSIZE_TOST, EFFSIZE_PWR.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'x1', @isnumeric);
addRequired(p, 'x2', @isnumeric);
addParameter(p, 'confLvl', 0.95, @(x) isnumeric(x) && x > 0 && x < 1);
addParameter(p, 'flgAbs', false, @islogical);
parse(p, x1, x2, varargin{:});

confLvl = p.Results.confLvl;
flgAbs  = p.Results.flgAbs;

%% ========================================================================
%  COMPUTE
%  ========================================================================

x1 = x1(~isnan(x1));
x2 = x2(~isnan(x2));
n1 = numel(x1);
n2 = numel(x2);

if n1 < 2 || n2 < 2
    error('effSize_d:tooFewSamples', ...
        'Each group needs at least 2 finite values (n1=%d, n2=%d).', n1, n2);
end

m1 = mean(x1); m2 = mean(x2);
s1 = std(x1);  s2 = std(x2);

sPooled = sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2));
d = (m1 - m2) / sPooled;

% Hedges' bias correction
J = 1 - 3 / (4 * (n1 + n2) - 9);
g = J * d;

% Standard error of d (Hedges & Olkin)
se_d = sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)));

% Confidence interval (normal approximation)
zCrit = norminv(0.5 + confLvl / 2);
ci = [d - zCrit * se_d, d + zCrit * se_d];

if flgAbs
    d = abs(d);
    g = abs(g);
    ci = sort(abs(ci));
end

%% ========================================================================
%  PACK INFO
%  ========================================================================

info = struct();
info.mean1   = m1;   info.mean2 = m2;
info.sd1     = s1;   info.sd2   = s2;
info.sPooled = sPooled;
info.n1      = n1;   info.n2    = n2;
info.se_d    = se_d;

end
