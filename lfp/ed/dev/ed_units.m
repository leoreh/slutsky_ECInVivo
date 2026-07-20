% ED_UNITS  Does the event produce the discharge's UNIT-FIRING signature?
%
% The one criterion here that is not circular. Every waveform measure
% describes the same trace the event was detected in, so it can only ever
% agree with the detector. Spiking is a physiological consequence, measured
% on different electrodes and a different signal.
%
% Maslarova et al. 2025 (Nat Commun), and their earlier chronic CA1 work:
% an IED drives a brief population burst and is then followed by a PROLONGED
% SUPPRESSION of firing, stronger in interneurons. A sharp wave-ripple drives
% a burst with no such suppression. So:
%
%   burst then long silence  -> discharge
%   burst, no silence        -> ripple / ordinary sharp wave
%   nothing                  -> artifact
%
% Classes are those of ed_laminar.m. Writes ed_units.mat + ed_units.png.

clear
DEVDIR = fileparts(mfilename('fullpath'));
BINSZ = 0.010;              % psth bin [s]
HALFW = 0.60;               % psth half-window [s]
BASE  = [-0.5 -0.2];        % baseline window [s]
SUPP  = [0.05 0.30];        % suppression window [s]
BURST = [-0.02 0.02];       % burst window [s]

load(fullfile(DEVDIR, 'ed_laminar.mat'), 'res');
bps = { ...
    'D:\Data\RA\raMCU1\raMCU1_080621_0930', ...
    'D:\Data\RA\raMCU2\raMCU2_080621_0930', ...
    'D:\Data\RA\raMCU3\raMCU3_211203_084720', ...
    'D:\Data\RA\raMCU4\raMCU4_211220_0834', ...
    'D:\Data\RA\raMCU5\raMCU5_220322_1906'};

edges = -HALFW : BINSZ : HALFW;
ctrs  = edges(1 : end - 1) + BINSZ / 2;
out = struct([]);

for iB = 1 : numel(bps)

    basepath = bps{iB};
    [~, basename] = fileparts(basepath);
    fprintf('\n===== %s =====\n', basename);

    S = load(fullfile(basepath, [basename '.spikes.cellinfo.mat']), 'spikes');
    st = S.spikes.times(:);
    nU = numel(st);
    fprintf('  %d units\n', nU);

    cls = {'ED', 'gate', 'rej', 'ripp', 'rand'};
    for iC = 1 : numel(cls)
        if ~isfield(res(iB), cls{iC}) || isempty(res(iB).(cls{iC}))
            continue
        end
        t = res(iB).(cls{iC}).t(:);
        if numel(t) < 5, continue, end

        % population psth: pooled over units, rate per unit per second
        cnt = zeros(numel(ctrs), nU);
        for iU = 1 : nU
            s = st{iU};
            for iE = 1 : numel(t)
                d = s(s > t(iE) - HALFW & s < t(iE) + HALFW) - t(iE);
                if ~isempty(d)
                    cnt(:, iU) = cnt(:, iU) + histcounts(d, edges)';
                end
            end
        end
        rate = cnt / (numel(t) * BINSZ);         % [nBin x nU] Hz

        iBase = ctrs >= BASE(1) & ctrs <= BASE(2);
        iSup  = ctrs >= SUPP(1) & ctrs <= SUPP(2);
        iBur  = ctrs >= BURST(1) & ctrs <= BURST(2);

        bl = mean(rate(iBase, :), 1);
        ok = bl > 0.05;                          % units with usable baseline
        gain = mean(rate(iBur, ok), 1) ./ bl(ok);
        supp = mean(rate(iSup, ok), 1) ./ bl(ok);

        pop = mean(rate(:, ok), 2);
        popBl = mean(pop(iBase));

        out(iB).(cls{iC}).pop  = pop / popBl;
        out(iB).(cls{iC}).gain = gain;
        out(iB).(cls{iC}).supp = supp;
        out(iB).(cls{iC}).n    = numel(t);

        fprintf(['  %-5s n=%3d  burst x%.2f   suppression x%.2f  ', ...
            '(%d/%d units suppressed)\n'], cls{iC}, numel(t), ...
            median(gain), median(supp), sum(supp < 0.8), sum(ok));
    end
    out(iB).basename = basename;
end

save(fullfile(DEVDIR, 'ed_units.mat'), 'out', 'ctrs');

% ------------------------------------------------------------------ figure
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
        plot(ctrs, c.pop, 'Color', clr.(cls{iC}), ...
            'LineWidth', 1.5, 'DisplayName', ...
            sprintf('%s (%d)', cls{iC}, c.n));
    end
    yline(1, ':'); xline(0, ':');
    title(out(iB).basename(1 : 6), 'Interpreter', 'none');
    xlabel('time (s)'); xlim([-0.4 0.5]);
    if iB == 1, ylabel('pop rate / baseline'); end
    legend('Location', 'northeast', 'Box', 'off', 'FontSize', 7);
end
title(tl, 'Peri-event population firing (discharge = burst then prolonged suppression)');
exportgraphics(fh, fullfile(DEVDIR, 'ed_units.png'), 'Resolution', 110);
fprintf('\nsaved ed_units.png\n');
