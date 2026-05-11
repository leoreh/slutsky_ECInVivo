function [out, info] = effSize_pwr(d, varargin)
% EFFSIZE_PWR Sample size or power for a two-sample independent t-test.
%
%   [out, info] = effSize_pwr(d, ...) computes either the per-group sample
%   size required to achieve a given power, or the achieved power for a
%   given per-group sample size. Mode is selected by which of 'n' or 'pwr'
%   is supplied; the other is solved for.
%
%   INPUT:
%       d           (double) Cohen's d (effect size; magnitude used).
%
%   OPTIONAL KEY-VALUE PAIRS (supply exactly one of n or pwr):
%       'n'         sample size  -> returns power. Scalar treated as
%                   equal n per group; 2-vector [n1, n2] uses the exact
%                   unequal-n formula.
%       'pwr'       desired power -> returns required n PER GROUP
%                   (assumes equal groups; returned as scalar).
%       'alpha'     two-sided significance level (default 0.05).
%       'flgExact'  use non-central t (true) or normal approximation
%                   (false). Default true. The two agree to within ~1 unit.
%
%   OUTPUT:
%       out         scalar. Achieved power (if n was given) or required
%                   per-group n (if pwr was given). n is rounded UP.
%       info        struct with d, alpha, n, pwr, mode ('n' | 'pwr'),
%                   and method ('exact' | 'approx').
%
%   Normal approximation:
%       n = 2 * (zAlpha + zBeta)^2 / d^2
%   Exact (non-central t):
%       power = 1 - nctcdf(tCrit, df, ncp) + nctcdf(-tCrit, df, ncp)
%       with ncp = d * sqrt(n/2) and df = 2*n - 2.
%
%   See also: EFFSIZE_D, EFFSIZE_TOST.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'd', @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'n', [], @(x) isempty(x) || (isnumeric(x) && all(x >= 2) && (isscalar(x) || numel(x) == 2)));
addParameter(p, 'pwr', [], @(x) isempty(x) || (isnumeric(x) && x > 0 && x < 1));
addParameter(p, 'alpha', 0.05, @(x) isnumeric(x) && x > 0 && x < 0.5);
addParameter(p, 'flgExact', true, @islogical);
parse(p, d, varargin{:});

n = p.Results.n;
pwr = p.Results.pwr;
alpha = p.Results.alpha;
flgExact = p.Results.flgExact;

if (isempty(n) && isempty(pwr)) || (~isempty(n) && ~isempty(pwr))
    error('effSize_pwr:badMode', ...
        'Provide exactly one of ''n'' or ''pwr'' (the other is solved for).');
end

dAbs = abs(d);
if dAbs == 0
    error('effSize_pwr:zeroEffect', 'd must be non-zero.');
end

%% ========================================================================
%  COMPUTE
%  ========================================================================

if ~isempty(n)
    % Mode: solve for power given n
    out = pwrFromN(n, dAbs, alpha, flgExact);
    mode = 'n';
    pwr = out;

else
    % Mode: solve for n given desired power
    if flgExact
        % Bisection on n (continuous, then ceil)
        lo = 2; hi = 1e5;
        while pwrFromN(hi, dAbs, alpha, true) < pwr && hi < 1e7
            hi = hi * 2;
        end
        for it = 1:60
            mid = 0.5 * (lo + hi);
            if pwrFromN(mid, dAbs, alpha, true) < pwr
                lo = mid;
            else
                hi = mid;
            end
            if hi - lo < 1e-3, break; end
        end
        out = ceil(hi);
    else
        zAlpha = norminv(1 - alpha / 2);
        zBeta  = norminv(pwr);
        out = ceil(2 * (zAlpha + zBeta)^2 / dAbs^2);
    end
    mode = 'pwr';
    n = out;
end

%% ========================================================================
%  PACK INFO
%  ========================================================================

info = struct();
info.d      = d;
info.alpha  = alpha;
info.n      = n;
info.pwr    = pwr;
info.mode   = mode;
if flgExact, info.method = 'exact'; else, info.method = 'approx'; end

end


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function pwr = pwrFromN(n, d, alpha, flgExact)

if isscalar(n), n1 = n; n2 = n; else, n1 = n(1); n2 = n(2); end
df = n1 + n2 - 2;
nEff = (n1 * n2) / (n1 + n2);       % harmonic-mean / 2; ncp scale factor

if flgExact
    ncp = d * sqrt(nEff);
    tCrit = tinv(1 - alpha / 2, df);
    pwr = 1 - nctcdf(tCrit, df, ncp) + nctcdf(-tCrit, df, ncp);
else
    zAlpha = norminv(1 - alpha / 2);
    pwr = 1 - normcdf(zAlpha - d * sqrt(nEff)) ...
            + normcdf(-zAlpha - d * sqrt(nEff));
end

end
