% ED_LAMINARPLOT  Per-channel mean ED waveform, with the detection channel
% marked. Reads ed_laminar.mat (see ed_laminar.m). Answers: does the discharge
% change shape across the probe, and where does the sleep_sig eeg average -
% the trace ED detection actually runs on - sit in that profile?

DEVDIR = fileparts(mfilename('fullpath'));
load(fullfile(DEVDIR, 'ed_laminar.mat'), 'res');

bps = { ...
    'D:\Data\RA\raMCU1\raMCU1_080621_0930', ...
    'D:\Data\RA\raMCU2\raMCU2_080621_0930', ...
    'D:\Data\RA\raMCU3\raMCU3_211203_084720', ...
    'D:\Data\RA\raMCU4\raMCU4_211220_0834', ...
    'D:\Data\RA\raMCU5\raMCU5_220322_1906'};

fh = figure('Position', [50 50 1800 900], 'Color', 'w');
tl = tiledlayout(1, numel(res), 'TileSpacing', 'compact');

for iB = 1 : numel(res)

    [~, bname] = fileparts(bps{iB});
    S = load(fullfile(bps{iB}, [bname '.sleep_sig.mat']), 'info');
    eegCh = S.info.eegCh;

    if isfield(res(iB), 'ED') && ~isempty(res(iB).ED)
        prof = res(iB).ED.prof;  lbl = 'curated ED';
    else
        prof = res(iB).gate.prof; lbl = 'gated';
    end
    ch = res(iB).chNeural;
    fs = res(iB).fs;
    tAx = ((1 : size(prof, 1)) - size(prof, 1) / 2) / fs * 1000;

    % one row per channel, spaced by a common offset
    off = 1.1 * max(abs(prof(:)));
    nexttile; hold on
    for iC = 1 : numel(ch)
        isDet = ismember(ch(iC), eegCh);
        if isDet, clr = [0.85 0.1 0.1]; lw = 2; else, clr = [0 0 0]; lw = 1; end
        plot(tAx, prof(:, iC) - (iC - 1) * off, 'Color', clr, ...
            'LineWidth', lw);
        text(tAx(1), -(iC - 1) * off, sprintf('%d ', ch(iC)), ...
            'HorizontalAlignment', 'right', 'FontSize', 8, 'Color', clr);
    end
    xline(0, ':');
    ampCh = max(abs(prof), [], 1);
    [~, iBest] = max(ampCh);
    title(sprintf('%s\n%s (n=%d)\ndetCh %s | best ch %d (%.0f uV)', ...
        bname(1 : 6), lbl, numel(res(iB).(ternary(isfield(res(iB), ...
        'ED') && ~isempty(res(iB).ED), 'ED', 'gate')).t), ...
        mat2str(eegCh), ch(iBest), ampCh(iBest)), ...
        'Interpreter', 'none', 'FontSize', 9);
    xlabel('time (ms)'); set(gca, 'YTick', []);
    xlim([-100 100]);

    fprintf('\n%s  detCh %s\n', bname, mat2str(eegCh));
    fprintf('  per-ch amp (uV): %s\n', ...
        strjoin(compose('%d:%.0f', ch(:), ampCh(:)'), '  '));
end

title(tl, ['Mean ED waveform per channel ', ...
    '(red = channels averaged for detection)']);
exportgraphics(fh, fullfile(DEVDIR, 'ed_laminar.png'), 'Resolution', 110);
fprintf('\nsaved ed_laminar.png\n');

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
