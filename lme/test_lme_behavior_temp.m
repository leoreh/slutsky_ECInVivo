
try
    tbl = table((1:10)', (1:10)' + randn(10,1), 'VariableNames', {'x', 'y'});
    fitlme(tbl, 'y ~ x');
    disp('fitlme worked without random effects');
catch ME
    disp(['fitlme failed: ' ME.message]);
end
