function y = stat_geomean(x)
% STAT_GEOMEAN Robust geometric mean for groupsummary
% Incorporates logic from tbl_tNorm for handling non-positive values

% Filter NaNs
x = x(~isnan(x));

if isempty(x)
    y = NaN;
    return;
end

% Handle non-positive values
posVals = x(x > 0);
if isempty(posVals)
    % If all values are <= 0, geomean is mathematically 0 (or undefined in log space)
    % Following tbl_tNorm logic, we might floor them, but if ALL are <=0, 
    % the "floor" based on positive values is undefined. 
    % Let's return 0 to be safe.
    y = 0;
    return;
end

% Determine floor value (half of min positive value)
floorVal = min(posVals) / 2;

% Apply floor to non-positive values
x(x <= 0) = floorVal;

% Calculate Geometric Mean
y = exp(mean(log(x)));

end
