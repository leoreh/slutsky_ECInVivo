function test_plot_raster()
% TEST_PLOT_RASTER Verifies functionality of graphics/plot_raster.m

fprintf('Testing plot_raster...\n');

% 1. Data Generation
nUnits = 5;
spktimes = cell(nUnits, 1);
for i = 1:nUnits
    % Random spikes
    spktimes{i} = sort(rand(1, 20 + i*5) * 10);
end

% 2. Test Different Plot Types

% --- Vertline ---
hFig1 = figure('Name', 'Test Vertline', 'Color', 'w');
[hAx1, hPlt1] = plot_raster(spktimes, 'hAx', gca, 'plotType', 'vertline', ...
    'color', 'k', 'vertSpikeH', 0.8, 'autoLabel', true);
title(hAx1, 'Vertline Test');

if isvalid(hPlt1)
    fprintf('  [PASS] Vertline plot created.\n');
else
    fprintf('  [FAIL] Vertline plot handle invalid.\n');
end

% --- Horzline ---
hFig2 = figure('Name', 'Test Horzline', 'Color', 'w');
[hAx2, hPlt2] = plot_raster(spktimes, 'hAx', gca, 'plotType', 'horzline', ...
    'color', 'b', 'spikeDur', 0.1, 'autoLabel', true);
title(hAx2, 'Horzline Test');

if isvalid(hPlt2)
    fprintf('  [PASS] Horzline plot created.\n');
else
    fprintf('  [FAIL] Horzline plot handle invalid.\n');
end

% --- Scatter ---
hFig3 = figure('Name', 'Test Scatter', 'Color', 'w');
[hAx3, hPlt3] = plot_raster(spktimes, 'hAx', gca, 'plotType', 'scatter', ...
    'color', 'r', 'markerSize', 10, 'autoLabel', true);
title(hAx3, 'Scatter Test');

if isvalid(hPlt3)
    fprintf('  [PASS] Scatter plot created.\n');
else
    fprintf('  [FAIL] Scatter plot handle invalid.\n');
end

fprintf('Please inspect figures 1, 2, and 3.\n');

end
