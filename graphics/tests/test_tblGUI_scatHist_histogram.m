function test_tblGUI_scatHist_histogram()
% Create dummy data
rng(42);
n = 1000;

% Log-normal data for X
x = 10.^(randn(n, 1) + 2);

% Normal data for Y
y = randn(n, 1) * 5 + 50;

% Groups
g = randi(3, n, 1);
gStr = "Group " + g;

tbl = table(x, y, gStr, 'VariableNames', {'ValX', 'ValY', 'Group'});

% Open GUI
hFig = tblGUI_scatHist(tbl, 'xVar', 'ValX', 'yVar', 'ValY', 'grpVar', 'Group');

% Access data to verify internals
data = hFig.UserData;

% 1. Test Linear Scale (Default)
% Trigger update (implicitly done by construction)

% Check Histogram properties
hHistX = findobj(data.hAxHistX, 'Type', 'histogram');
if isempty(hHistX)
    error('Histogram X not found');
end

disp('Checking Linear Scale Histogram Properties...');
check_hist_props(hHistX, 'linear');

% Check for Lines (KDE)
hLinesX = findobj(data.hAxHistX, 'Type', 'line');
if isempty(hLinesX)
    error('KDE lines not found in X axis');
end
fprintf('Found %d KDE lines in X axis (Linear)\n', length(hLinesX));


% 2. Change X to Log Scale
fprintf('Switching X to Log Scale...\n');
set(data.ddXScale, 'Value', 2); % 2 = Log
% Dispatch callback
cbk = get(data.ddXScale, 'Callback');
cbk(data.ddXScale, []);

% Check again
hHistX = findobj(data.hAxHistX, 'Type', 'histogram');
disp('Checking Log Scale Histogram Properties...');
check_hist_props(hHistX, 'log');

hLinesX = findobj(data.hAxHistX, 'Type', 'line');
if isempty(hLinesX)
    error('KDE lines not found in X axis (Log)');
end
fprintf('Found %d KDE lines in X axis (Log)\n', length(hLinesX));


% 3. Check Y Histogram (Linear)
hHistY = findobj(data.hAxHistY, 'Type', 'histogram');
disp('Checking Y Scale Histogram Properties (Linear)...');
check_hist_props(hHistY, 'linear');

hLinesY = findobj(data.hAxHistY, 'Type', 'line');
if isempty(hLinesY)
    error('KDE lines not found in Y axis');
end

% 4. Test Discretize Feature
disp('Testing Discretize Feature...');
% Enable Discretize
set(data.chkDisc, 'Value', 1);
% Trigger update
cbk = get(data.chkDisc, 'Callback');
cbk(data.chkDisc, []);

% Check for ErrorBars
hErr = findall(data.hAxScatter, 'Type', 'errorbar');
if isempty(hErr)
    error('ErrorBars not found after enabling Discretize');
end
fprintf('Found %d ErrorBar objects.\n', length(hErr));

% Check Scatter Transparency
hScat = findobj(data.hAxScatter, 'Type', 'scatter');
if isempty(hScat), error('Scatter plot missing'); end

% Scatter has Multiple objects (one per group), check first one
alphaVal = hScat(1).MarkerFaceAlpha;
% Original default is in data.defAlpha (0.6). We expect 0.6 * 0.3 = 0.18
expectedAlpha = data.defAlpha * 0.3;
if abs(alphaVal - expectedAlpha) > 0.05
    error('Scatter Alpha %.2f does not match expected %.2f', alphaVal, expectedAlpha);
end

% 5. Test Deactivation (Cleanup)
disp('Testing Deactivation (Cleanup)...');
set(data.chkDisc, 'Value', 0);
cbk(data.chkDisc, []);

hErr = findall(data.hAxScatter, 'Type', 'errorbar');
if ~isempty(hErr)
    error('ErrorBars PERSISTED after disabling Discretize! Found %d.', length(hErr));
else
    disp('ErrorBars successfully removed.');
end

fprintf('Verification Passed!\n');
close(hFig);
end

function check_hist_props(hHistArray, scaleType)
% hHistArray might be multiple histograms (one per group)
for i = 1:length(hHistArray)
    h = hHistArray(i);

    % Check Style - Should now be BAR not stairs
    if ~strcmpi(h.DisplayStyle, 'bar')
        error('DisplayStyle is not bar. Got: %s', h.DisplayStyle);
    end

    % Check Normalization - Should now be PDF
    if ~strcmpi(h.Normalization, 'pdf')
        error('Normalization is not pdf. Got: %s', h.Normalization);
    end

    % Check FaceAlpha
    if h.FaceAlpha > 0.6
        warning('FaceAlpha might be too high: %.2f', h.FaceAlpha);
    end

    % Check Bin Edges
    edges = h.BinEdges;
    diffs = diff(edges);

    if strcmp(scaleType, 'linear')
        % Diffs should be constant (approx)
        if std(diffs) > 1e-6
            error('Bin edges are not linear for linear scale.');
        end
    elseif strcmp(scaleType, 'log')
        % Log of edges should be linear -> diff of log edges constant
        logEdges = log10(edges);
        diffLog = diff(logEdges);
        if std(diffLog) > 1e-6
            error('Bin edges are not log-spaced for log scale.');
        end
    end
end
end
