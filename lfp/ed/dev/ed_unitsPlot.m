% ED_UNITSPLOT  Draw the peri-event population firing from ed_units.mat.

DEVDIR = fileparts(mfilename('fullpath'));
load(fullfile(DEVDIR, 'ed_units.mat'), 'out', 'ctrs');

fh = figure('Position', [50 50 1700 400], 'Color', 'w');
tl = tiledlayout(1, numel(out), 'TileSpacing', 'compact');
clr = struct('ED', [0.8 0 0], 'gate', [0.9 0.5 0], 'rej', [0.4 0.4 0.4], ...
    'ripp', [0 0.4 0.8], 'rand', [0.7 0.7 0.7]);

for iB = 1 : numel(out)
    nexttile; hold on
    cls = intersect(fieldnames(out(iB)), fieldnames(clr), 'stable');
    for iC = 1 : numel(cls)
        c = out(iB).(cls{iC});
        if isempty(c), continue, end
        plot(ctrs, c.pop, 'Color', clr.(cls{iC}), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%s (%d)', cls{iC}, c.n));
    end
    yline(1, ':'); xline(0, ':');
    title(out(iB).basename(1 : 6), 'Interpreter', 'none');
    xlabel('time (s)'); xlim([-0.4 0.5]); ylim([0 2]);
    if iB == 1, ylabel('pop rate / baseline'); end
    legend('Location', 'southeast', 'Box', 'off', 'FontSize', 7);
end
title(tl, ['Peri-event population firing  ', ...
    '(discharge = prolonged suppression; ripple / artifact = none)']);
exportgraphics(fh, fullfile(DEVDIR, 'ed_units.png'), 'Resolution', 110);
fprintf('saved ed_units.png\n');
