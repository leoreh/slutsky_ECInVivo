function test_plot_axSize_mod()
% TEST_PLOT_AXSIZE_MOD Verifies flgFullscreen and flgPos logic.

fprintf('Testing plot_axSize...\n');

% 1. Test flgPos (offset)
fprintf('  Launching Figure 1 (flgPos)... check if offset on screen.\n');
hFig1 = plot_axSize('flgPos', true, 'axShape', 'square', 'szOnly', false);
title(gca, 'flgPos = true');
pause(1);

% 2. Test flgFullscreen
fprintf('  Launching Figure 2 (flgFullscreen)... check if fullscreen on 2nd monitor (if exists).\n');
hFig2 = plot_axSize('flgFullscreen', true, 'axShape', 'wide', 'szOnly', false);
title(gca, 'flgFullscreen = true');

% Print monitor info
monitors = get(0, 'MonitorPositions');
fprintf('  Detected %d monitor(s).\n', size(monitors, 1));
fprintf('  Fig 1 Position: %s\n', mat2str(hFig1.Position));
fprintf('  Fig 2 Position: %s\n', mat2str(hFig2.Position));

fprintf('Test complete. Close figures manually.\n');

end
