function d = norm_xcorr_dist2(x,y)
    
    arguments (Input)
        x (:,1) {mustBeNumeric}
        y (:,1) {mustBeNumeric}
    end
    arguments (Output)
        d (1,1) {mustBeNumeric}
    end
    
    x = x-mean(x);
    y = y-mean(y);
    [r,lags] = xcorr(x,y,"normalized");
    
    % taken from
    % https://rdrr.io/cran/TSdist/src/R/cross_correlation_distance.R.
    % note that only sum over lags<0, I dont know why...
    d = sqrt( ...
        (1 - round(r(lags==0).^2,5)) ./ ...
        sum(r(lags<0).^2) ...
        );
end