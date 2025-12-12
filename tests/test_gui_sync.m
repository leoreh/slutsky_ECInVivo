function test_gui_sync()
% TEST_GUI_SYNC Reproduces synchronization issues in utypes_gui components.
%
% Steps:
% 1. Creates a dummy table.
% 2. Opens utypes_gui (mocked).
% 3. Programmatically changes Grouping in one window.
% 4. Checks if other windows updated.
% 5. Mocks a selection and checks callbacks.

clc; close all;

%% 1. Setup Data
nPoints = 100;
tAxis = linspace(0, 1, 50);
uTbl = table();
uTbl.Name = repmat({'Test'}, nPoints, 1);
uTbl.UnitID = (1:nPoints)';
uTbl.Grp = categorical(randi([1 3], nPoints, 1));
uTbl.XVal = randn(nPoints, 1);
uTbl.YVal = randn(nPoints, 1) + 5;
uTbl.UnitType = uTbl.Grp; % Secondary group var

% Create generic traces
uTbl.Trace = num2cell(randn(nPoints, 50), 2);

% Cfg
cfg = struct();
cfg.xVar = 'XVal';
cfg.yVar = 'YVal';
cfg.grpVar = 'Grp';

fprintf('Launching utypes_gui...\n');

% Mock basepaths empty as we provide uTbl
hFigMain = utypes_gui('uTbl', uTbl, 'cfg', cfg, 'tAxis', tAxis);
pause(1); % Allow GUI to render

% Get Handles of all figures
figs = findall(0, 'Type', 'figure');
hScat = []; hTrace = [];
for i=1:length(figs)
    if contains(figs(i).Name, 'Scatter Plot'), hScat = figs(i); end
    if contains(figs(i).Name, 'Traces'), hTrace = figs(i); end
end

if isempty(hScat) || isempty(hTrace)
    error('Failed to find GUI figures.');
end

%% 2. Test Group Synchronization (Scatter -> Traces)
fprintf('\n--- TEST 1: Group Sync (Scatter -> Traces) ---\n');

% Initial State
dScat = hScat.UserData;
dTrace = hTrace.UserData;
fprintf('Initial Scat Group: %s\n', get_selected_string(dScat.ddGrp));
fprintf('Initial Trace Group: %s\n', get_selected_string(dTrace.ddGrpBy));

% Change Group in Scatter Programmatically
fprintf('Changing Scatter Group to "UnitType"...\n');
set_dropdown(dScat.ddGrp, 'UnitType');
% Trigger Callback manually (as if user selected)
dScat.ddGrp.Callback(dScat.ddGrp, []);

pause(0.5);

% Check Trace Window
dTrace = hTrace.UserData; % Reload
currTraceGrp = get_selected_string(dTrace.ddGrpBy);
fprintf('Current Trace Group: %s\n', currTraceGrp);

if strcmp(currTraceGrp, 'UnitType')
    fprintf('PASSED: Trace window updated to "UnitType".\n');
else
    fprintf('FAILED: Trace window is "%s", expected "UnitType".\n', currTraceGrp);
end

%% 3. Test Selection Sync (Selection Box)
fprintf('\n--- TEST 2: Selection Box Update ---\n');
% We can't easily draw a polygon programmatically without lower level tools,
% but we can call the selection callback directly.
% Actually, the issue reported was "selecting specific boxes doesn't change the scatter or historgam plots".
% This likely refers to the Checkboxes in the side panel.

% Let's test Checkbox filtering.
dScat = hScat.UserData;
chk = dScat.chkGrp;
if isempty(chk)
    fprintf('FAILED: No checkboxes found in Scatter GUI.\n');
else
    fprintf('Found %d checkboxes.\n', length(chk));
    % Uncheck first one
    fprintf('Unchecking category: %s\n', chk(1).String);
    chk(1).Value = 0;
    % Trigger Callback
    chk(1).Callback(chk(1), []);

    pause(0.5);

    % Verify Plot Data
    % We need to check if points were removed.
    % In tblGUI_scatHist, scatter points have HandleVisibility 'off' if I recall correctly from reading...
    % Wait, let's check code.
    % line 482: 'HitTest', 'off', 'PickableParts', 'none'
    % But they are children of hAxScatter.

    kids = dScat.hAxScatter.Children;
    % Expected: Number of scatter objects should equal number of ACTIVE groups.
    % If we unchecked one, we expect one less group displayed (if it filters groups entirely)
    % OR the same number of groups but one has 0 points (or is not plotted).
    % The code loop iterates 1:length(grpLabels).
    % Inside the loop: `if isGrpActive && ~isempty(data.chkGrp)` -> filters grpLabels.

    % So we expect fewer scatter objects?
    numScatter = 0;
    for k=1:length(kids)
        if isa(kids(k), 'matlab.graphics.chart.primitive.Scatter')
            numScatter = numScatter + 1;
        end
    end

    fprintf('Number of Scatter objects found: %d\n', numScatter);

    % We have 3 groups (randi([1 3])). Unchecking 1 should leave 2.
    if numScatter == 2
        fprintf('PASSED: Scatter plot updated (2 groups visible).\n');
    else
        fprintf('FAILED: Expected 2 scatter groups, found %d.\n', numScatter);
    end
end

fprintf('\nDone.\n');

end

function str = get_selected_string(hObj)
items = get(hObj, 'String');
val = get(hObj, 'Value');
str = items{val};
end

function set_dropdown(hObj, str)
items = get(hObj, 'String');
idx = find(strcmp(items, str));
if ~isempty(idx)
    set(hObj, 'Value', idx);
else
    warning('Item not found: %s', str);
end
end
