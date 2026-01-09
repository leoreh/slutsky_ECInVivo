function test_lme_ablation_comb()

% Create dummy data: 10 subjects, 5 per group
nSub = 10;
nObsPerSub = 10; % 10 observations per subject

names = {};
groups = {};
y = [];
x = [];

% Group A: Sub1-5
for i = 1:5
    nm = sprintf('SubA_%d', i);
    names = [names; repmat({nm}, nObsPerSub, 1)];
    groups = [groups; repmat({'A'}, nObsPerSub, 1)];
    y = [y; randn(nObsPerSub, 1)];
    x = [x; randn(nObsPerSub, 1)];
end

% Group B: Sub6-10
for i = 1:5
    nm = sprintf('SubB_%d', i);
    names = [names; repmat({nm}, nObsPerSub, 1)];
    groups = [groups; repmat({'B'}, nObsPerSub, 1)];
    y = [y; randn(nObsPerSub, 1)];
    x = [x; randn(nObsPerSub, 1)];
end

tbl = table(names, groups, x, y, 'VariableNames', {'Name', 'Group', 'X', 'Y'});

% Formula: Y ~ X + Group + (1|Name)
frml = 'Y ~ X + Group + (1|Name)';

fprintf('Running lme_ablation with 5x5 groups...\n');

% We expect 25 iterations regardless of nReps (if exhaustive works).
% If it fails, it will use nReps=2 -> 10 iterations.
res = lme_ablation(tbl, frml, 'dist', 'Normal', 'nReps', 2, 'flgPlot', false);

nIterActual = size(res.perfPrim, 1);
fprintf('Actual Iterations: %d\n', nIterActual);

if nIterActual == 25
    fprintf('SUCCESS: 25 iterations found as expected for 5x5 design.\n');
else
    fprintf('FAILURE: Expected 25 iterations, got %d.\n', nIterActual);
end

end
