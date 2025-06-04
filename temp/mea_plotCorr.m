% Helper function for plotting correlations
function mea_plotCorr(x, y, flg_log)

% Remove outliers from both x and y
% [~, TFx] = rmoutliers(x(:));
% [~, TFy] = rmoutliers(y(:));
% % Remove points where either x or y was an outlier
% valid = ~(TFx | TFy);
% x = x(valid);
% y = y(valid);

x = x(:);
y = y(:);

scatter(x, y, 'filled');
hold on;

if flg_log
    set(gca, 'xscale', 'log');
    % set(gca, 'yscale', 'log');
    valid_idx = x > 0 & ~isnan(x) & ~isnan(y) & ~isinf(y) & ~isinf(x);
    if sum(valid_idx) > 1
        r = corrcoef(log10(x(valid_idx)), y(valid_idx));
        text(0.05, 0.95, sprintf('r = %.3f (n=%d)', r(1,2), sum(valid_idx)), 'Units', 'normalized');
    end
else
    valid_idx = ~isnan(x) & ~isnan(y);
    if sum(valid_idx) > 1
        r = corrcoef(x(valid_idx), y(valid_idx));
        text(0.05, 0.95, sprintf('r = %.3f (n=%d)', r(1,2), sum(valid_idx)), 'Units', 'normalized');
    end
end

grid on;
end

