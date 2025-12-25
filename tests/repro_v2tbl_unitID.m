% Reproduction script for non-deterministic unitID in v2tbl
% Scenario 1: Default behavior (Persistent variable) -> IDs should increment
% Scenario 2: Explicit offset -> IDs should be deterministic

% Create dummy data
v(1).data = [1; 2; 3];
v(2).data = [4; 5];

% Create varMap
varMap.Value = 'data';

%% Test 1: Default Behavior (Persistence)
disp('--- Test 1: Default Behavior (Persistence) ---');
% Call v2tbl first time (assuming fresh state or continuation)
tbl1 = v2tbl('v', v, 'varMap', varMap);
firstIDs = tbl1.UnitID';
fprintf('Run 1 UnitIDs: %s\n', mat2str(firstIDs));

% Call v2tbl second time - should continue incrementing
tbl2 = v2tbl('v', v, 'varMap', varMap);
secondIDs = tbl2.UnitID';
fprintf('Run 2 UnitIDs: %s\n', mat2str(secondIDs));

% Check that they are DIFFERENT (and specifically, second > first)
if all(secondIDs > max(firstIDs))
    disp('SUCCESS: Default behavior preserves varying UnitIDs.');
else
    disp('FAILURE: Default behavior did not increment UnitIDs as expected.');
    % If we haven't fixed it yet, this might actually return SUCCESS effectively if the persistent var works
    % But if we just cleared functions or something, it might fail.
    % For the purpose of "repro script", the current codebase DOES increment.
end

%% Test 2: Deterministic Behavior (Explicit Offset)
disp('--- Test 2: Deterministic Behavior (uOffset = 0) ---');
% Call with explicit offset 0
tbl3 = v2tbl('v', v, 'varMap', varMap, 'uOffset', 0);
thirdIDs = tbl3.UnitID';
fprintf('Run 3 (Offset 0) UnitIDs: %s\n', mat2str(thirdIDs));

% Call again with explicit offset 0
tbl4 = v2tbl('v', v, 'varMap', varMap, 'uOffset', 0);
fourthIDs = tbl4.UnitID';
fprintf('Run 4 (Offset 0) UnitIDs: %s\n', mat2str(fourthIDs));

% Check that they are EQUAL
if isequal(tbl3.UnitID, tbl4.UnitID)
    disp('SUCCESS: Explicit uOffset produces deterministic UnitIDs.');
else
    disp('FAILURE: Explicit uOffset did NOT produce deterministic UnitIDs.');
end

% Bonus: Check if Run 3 matches what we expect for a "fresh" run (UnitIDs starting roughly at 1000 + row)
% v2tbl logic: uOffset + (iFile * 1000) + rowIdx
% File 1: 0 + 1000 + [1,2,3] -> [1001, 1002, 1003]
% File 2: 0 + 2000 + [1,2] -> [2001, 2002]
expectedIDs = [1001, 1002, 1003, 2001, 2002]';
if isequal(tbl3.UnitID, expectedIDs)
    disp('SUCCESS: UnitIDs match expected calculation logic.');
else
    disp('FAILURE: UnitIDs do not match expected logic.');
    fprintf('Expected: %s\n', mat2str(expectedIDs'));
    fprintf('Actual:   %s\n', mat2str(thirdIDs));
end
