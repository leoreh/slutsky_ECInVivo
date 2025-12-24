% Test script for tblGUI_raster

%% Create Dummy Table
nUnits = 200; % Increased
spktimes = cell(nUnits, 1);
group = strings(nUnits, 1);

for i = 1:nUnits
    % Random spikes up to 4 hours (14400s)
    spktimes{i} = sort(rand(1, 100) * 14400);
    if i <= 90
        group(i) = "Ctrl1";
    elseif i <= 180
        group(i) = "Ctrl2";
    else
        group(i) = "Treat";
    end
end

tbl = table(spktimes, categorical(group), 'VariableNames', {'spktimes', 'Group'});

%% Run GUI
close all;
fprintf('Running tblGUI_raster with grpVal="Ctrl2"...\n');

try
    hFig = tblGUI_raster(tbl, 'grpVar', 'Group', 'grpVal', 'Ctrl2');
    title(hFig.Children(end), 'Test Raster GUI (Ctrl2 Only)');
    fprintf('GUI launched successfully.\n');

    % Verify filtering visually:
    % Should see nUnits = 90 (rows 91-180), grouped compactly 1-90.
    % If logic is wrong, it might show 91-180.

catch ME
    fprintf('Error launching GUI: %s\n', ME.message);
    rethrow(ME);
end

%% Cleanup
% close(hFig);
