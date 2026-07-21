% ED_ALIGNDIAG  Where does t = 0 actually sit on the raw deflection?
%
% Leore's cluster averages show a NOTCH at t = 0 - a local maximum - with the
% real trough 5-7 ms later, in every cluster. That is the signature of aligning
% on a band-passed peak: 60-150 Hz rings with a 7-17 ms period, so max|filt|
% can land half a cycle off the extremum, and aligning every event to the same
% lobe imprints the ringing on the average.
%
% Measures, per mouse:
%   1. offset from t = 0 to the raw extremum, and whether it depends on polarity
%   2. how much of the pool lands further than half a ring period away
%   3. what the cluster-average waveform looks like before and after refining
%   4. whether the scalar features re-introduce the amplitude that
%      normalisation removes (4 of the 5 are monotone in size)

clear
DEVDIR = fileparts(mfilename('fullpath'));
MAPDUR = [-0.1 0.1];
SRCH = 0.012;                       % refinement search half-window [s]

bps = { ...
    'D:\Data\RA\raMCU3\raMCU3_211203_084720', ...
    'D:\Data\RA\raMCU4\raMCU4_211220_0834', ...
    'D:\Data\RA\raMCU5\raMCU5_220322_1906'};

met = ed_methods('default');

fprintf('\n%-8s %6s | %7s %7s %7s | %7s %7s\n', 'mouse', 'pool', ...
    'med dt', 'p25', 'p75', '|dt|>4ms', 'pol dif');
D = struct([]);
for iB = 1 : numel(bps)
    [ed, aux] = ed_detect(bps{iB}, 'met', met);
    [~, bname] = fileparts(bps{iB});
    maps = evt_maps(struct('lfp', aux.edSig.lfp, 'filt', aux.edSig.filt), ...
        ed.peakTime, aux.fs, 'mapDur', MAPDUR);

    keep = evt_gate(ed, met.qa);
    wv = double(maps.lfp(keep, :));
    ft = double(maps.filt(keep, :));
    tst = maps.tstamps;

    [dt, pol] = peakShift(wv, tst, SRCH);
    fprintf('%-8s %6d | %7.1f %7.1f %7.1f | %6.0f%% %7.1f\n', ...
        bname(1 : 6), nnz(keep), median(dt) * 1000, ...
        prctile(dt, 25) * 1000, prctile(dt, 75) * 1000, ...
        mean(abs(dt) > 0.004) * 100, ...
        (median(dt(pol > 0)) - median(dt(pol < 0))) * 1000);

    D(iB).bname = bname;
    D(iB).wv = wv; D(iB).ft = ft; D(iB).tst = tst; D(iB).dt = dt;
    D(iB).sc = [ed.fastZ(keep), ed.isoZ(keep), ed.posZ(keep), ...
        ed.amp(keep), ed.dur(keep)];
end

%% ========================================================================
%  ARE THE SCALARS JUST AMPLITUDE?
%  ========================================================================
% Normalisation strips size from the waveform; the scalar block can hand it
% straight back. Spearman against |amp|, which is size by definition.
nm = {'fastZ', 'isoZ', 'posZ', 'amp', 'dur'};
fprintf('\n%-8s |', 'rho vs |amp|');
fprintf(' %7s', nm{:}); fprintf('\n');
for iB = 1 : numel(D)
    r = corr(abs(D(iB).sc(:, 4)), D(iB).sc, 'type', 'Spearman', ...
        'rows', 'pairwise');
    fprintf('%-8s |', D(iB).bname(1 : 6));
    fprintf(' %7.2f', r); fprintf('\n');
end

%% ========================================================================
%  BEFORE / AFTER
%  ========================================================================
hFig = figure('Position', [80 80 1400 700]);
tl = tiledlayout(hFig, 2, numel(D), 'TileSpacing', 'compact');
for iPass = 1 : 2
    for iB = 1 : numel(D)
        if iPass == 1
            wv = D(iB).wv; tst = D(iB).tst; ttl = 'detected (max|filt|)';
        else
            [wv, tst] = realign(D(iB).wv, D(iB).tst, D(iB).dt);
            ttl = 'refined (raw extremum)';
        end
        cid = evt_clust(wv, tst, 'win', met.clust.win, 'nPC', met.clust.nPC, ...
            'nClust', met.clust.nClust, 'scalar', D(iB).sc, ...
            'detrend', met.clust.detrend, 'norm', met.clust.norm);

        ax = nexttile(tl); hold(ax, 'on');
        cnt = accumarray(cid(~isnan(cid)), 1);
        [~, ord] = sort(cnt, 'descend');
        for iK = 1 : min(3, numel(ord))
            m = cid == ord(iK);
            plot(ax, tst * 1000, median(wv(m, :), 1), 'LineWidth', 2, ...
                'DisplayName', sprintf('%d (n=%d)', ord(iK), nnz(m)));
        end
        xline(ax, 0, 'k:');
        xlim(ax, [-60 60]); grid(ax, 'on');
        title(ax, sprintf('%s - %s', D(iB).bname(1 : 6), ttl));
        legend(ax, 'Location', 'southwest');
    end
end
xlabel(tl, 'time (ms)'); ylabel(tl, 'lfp');
exportgraphics(hFig, fullfile(DEVDIR, 'ed_alignDiag.png'), 'Resolution', 150);
fprintf('\nsaved ed_alignDiag.png\n');


% =========================================================================
%  LOCAL
% =========================================================================
function [dt, pol] = peakShift(wv, tst, srch)
% Offset from t = 0 to the raw extremum within SRCH, after a flank-fitted line
% is removed - otherwise a slow wander sets the extremum, not the event.
t = tst(:);
iEdge = abs(t) >= 0.5 * max(abs(t));
A = [t(iEdge) - mean(t(iEdge)), ones(nnz(iEdge), 1)];
X = wv - ([t - mean(t(iEdge)), ones(numel(t), 1)] * (A \ wv(:, iEdge)'))';

iS = abs(t) <= srch;
[~, iPk] = max(abs(X(:, iS)), [], 2);
tS = t(iS);
dt = tS(iPk);
i0 = find(iS, 1);
pol = sign(X(sub2ind(size(X), (1 : size(X, 1))', i0 + iPk - 1)));

end     % peakShift


function [wv, tst] = realign(wv, tst, dt)
% Re-cut each waveform so t = 0 sits at its own raw extremum. The map is wider
% than the clustering window, so this is an index shift, not a re-read.
fs = 1 / median(diff(tst));
sh = round(dt(:) * fs);
keep = abs(tst) <= max(abs(tst)) - max(abs(sh)) / fs;
i0 = find(keep);
out = nan(size(wv, 1), numel(i0));
for iEv = 1 : size(wv, 1)
    idx = i0 + sh(iEv);
    ok = idx >= 1 & idx <= size(wv, 2);
    out(iEv, ok) = wv(iEv, idx(ok));
end
wv = out;
tst = tst(i0);
tst = tst - tst(round(end / 2));

end     % realign
