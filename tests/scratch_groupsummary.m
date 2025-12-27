
% Synthetic Data
tbl = table();
tbl.Group = {'A'; 'A'; 'B'; 'B'; 'A'};
tbl.Data = [1 2; 3 4; 5 6; 7 8; 9 10]; % A: [1 2, 3 4, 9 10], B: [5 6, 7 8]

% Expected for A (pooled): mean([1 2 3 4 9 10]) = 4.8333
% Expected for B (pooled): mean([5 6 7 8]) = 6.5

% Try groupsummary with custom function
disp('--- Custom Function with Pooling ---')
try
    % Note: x will be a matrix of size [RowsInGroup, Cols]
    statFunc = @(x) mean(x, 'all');
    S = groupsummary(tbl, 'Group', statFunc, 'Data')
catch ME
    disp(ME.message)
end

% Try groupsummary with standard mean (check column-wise)
disp('--- Standard Mean (Column-wise?) ---')
try
    S2 = groupsummary(tbl, 'Group', 'mean', 'Data')
    % If column wise, A should be [mean([1 3 9]), mean([2 4 10])] = [4.33, 5.33]
catch ME
    disp(ME.message)
end

% Try multiple outputs in custom function
disp('--- Custom Function Multiple Outputs ---')
try
    % Returns [mean, std]
    statFunc2 = @(x) [mean(x, 'all'), std(x, 0, 'all')];
    S3 = groupsummary(tbl, 'Group', statFunc2, 'Data')
catch ME
    disp(ME.message)
end
