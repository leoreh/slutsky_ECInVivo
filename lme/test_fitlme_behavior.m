
% Create dummy data
tbl = table(randn(100,1), randn(100,1), randn(100,1), 'VariableNames', {'Y', 'A', 'B'});

try
    % Try fitting with NO random effects using fitlme
    mdl = fitlme(tbl, 'Y ~ A + B');
    disp('SUCCESS: fitlme accepted formula without random effects.');
    disp(['Class of mdl: ' class(mdl)]);

    % Check if fixedEffects works
    [fe, names, stats] = fixedEffects(mdl);
    disp('SUCCESS: fixedEffects method works on the result.');
    disp(stats);

catch ME
    disp('FAILURE: fitlme threw an error.');
    disp(ME.message);
end
