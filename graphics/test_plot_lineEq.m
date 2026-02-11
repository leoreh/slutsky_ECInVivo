
try
    figure('Color', 'w');
    subplot(2,2,1);
    scatter(1:10, 1:10);
    plot_lineEq();
    title('Default');

    subplot(2,2,2);
    scatter(randn(100,1), randn(100,1));
    h = plot_lineEq('lineStyle', '-r', 'flgZero', false);
    title('Random Data, Red Line, No Zero');

    subplot(2,2,3);
    scatter([1, 100], [1, 50]);
    plot_lineEq('lims', 150);
    title('Explicit Lims');

    subplot(2,2,4);
    plot_lineEq();
    title('Empty Axes');
    
    fprintf('Test Passed!\n');
catch ME
    fprintf('Test Failed: %s\n', ME.message);
    rethrow(ME);
end
