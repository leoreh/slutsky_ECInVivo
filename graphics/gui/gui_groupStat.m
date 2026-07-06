function [m, lo, hi] = gui_groupStat(vals, statType, varargin)
% GUI_GROUPSTAT  Central tendency and error bounds, computed down dim 1.
%
%   [m, lo, hi] = gui_groupStat(vals, statType) treats the rows of VALS as
%   observations and returns the center M and the absolute lower / upper
%   bounds LO / HI (1 x size(vals,2)). VALS may be a column vector (scalar
%   per group, as in guiTbl_bar) or a matrix (one value per x-sample, as in
%   guiTbl_xy).
%
%   statType (case-insensitive):
%       'Arithmetic'  mean +/- SEM
%       'Geometric'   geometric mean, multiplicative SEM
%       'Median'      median +/- notch (1.57*IQR/sqrt(n))
%
%   Name-value:
%       'Floor'  for Geometric, clamp values up to this floor instead of
%                dropping non-positive ones. {[]} (drop, matching guiTbl_bar)
%
%   Callers convert to whatever they need: error-bar half-widths are
%   (m - lo) and (hi - m); a shaded band uses lo and hi directly.

p = inputParser;
addParameter(p, 'Floor', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
parse(p, varargin{:});
floorVal = p.Results.Floor;

switch lower(statType)
    case 'geometric'
        v = vals;
        if isempty(floorVal)
            v(v <= 0) = NaN;            % drop non-positive (omitted below)
        else
            v = max(v, floorVal);       % clamp to detection floor
        end
        logv = log(v);
        n    = sum(~isnan(logv), 1);
        mLog = mean(logv, 1, 'omitnan');
        sLog = std(logv, 0, 1, 'omitnan') ./ sqrt(n);
        m  = exp(mLog);
        lo = exp(mLog - sLog);
        hi = exp(mLog + sLog);

    case 'median'
        m     = median(vals, 1, 'omitnan');
        q1    = prctile(vals, 25, 1);
        q3    = prctile(vals, 75, 1);
        n     = sum(~isnan(vals), 1);
        notch = 1.57 * (q3 - q1) ./ sqrt(n);
        lo = m - notch;
        hi = m + notch;

    otherwise   % arithmetic
        n = sum(~isnan(vals), 1);
        m = mean(vals, 1, 'omitnan');
        s = std(vals, 0, 1, 'omitnan') ./ sqrt(n);
        lo = m - s;
        hi = m + s;
end

end     % EOF
