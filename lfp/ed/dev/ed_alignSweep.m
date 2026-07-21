% ED_ALIGNSWEEP  Does baseline, alignment or normalisation limit the clusters?
%
% Three questions, all of them about what the waveform looks like BEFORE the
% components are taken:
%
%   DETREND. ed_clust fits its linear baseline over the whole window, the event
%     included, so a large asymmetric deflection tilts the line it is measured
%     against - and tilts it by an amount that depends on the event's own
%     polarity and asymmetry. snipFromBinary fits on the flanks instead. Does
%     it matter here?
%
%   ALIGNMENT. Detection puts t = 0 at max|filt|, the peak of the 60-150 Hz
%     trace. That is the sharpest point, not the raw extremum, and the offset
%     between them need not be the same for an upward and a downward event. If
%     it is not, the pool is misaligned by polarity and no amount of clustering
%     fixes it.
%
%   NORMALISATION. Unit peak (L-inf, current) versus L2 versus none.
%
% Also answers a prior question the sweep cannot: are the curated discharges of
% ONE mouse all the same polarity? If they are, polarity is a property of the
% channel and flipping events would merge two real types. If they are not,
% polarity is incidental and a sign-blind representation is worth testing.
%
% Scored as every other sweep here: events you must review to reach 80% of the
% curated discharges, taking clusters in order of discharge density.

clear
DEVDIR = fileparts(mfilename('fullpath'));
TOL = 0.050;
MAPDUR = [-0.1 0.1];
SRCH = 0.010;                       % realignment search half-window [s]

bps = { ...
    'D:\Data\RA\raMCU3\raMCU3_211203_084720', ...
    'D:\Data\RA\raMCU4\raMCU4_211220_0834', ...
    'D:\Data\RA\raMCU5\raMCU5_220322_1906'};

DT   = {'full', 'edge', 'none'};
NM   = {'peak', 'l2', 'none'};
ALGN = {'filt', 'raw'};
KS   = [20 30 45];

met = ed_methods('default');
cur = load(fullfile(DEVDIR, 'ed_curatedTimes.mat'));
cur = cur.edCurated;

%% ========================================================================
%  LOAD + POLARITY
%  ========================================================================
D = struct([]);
fprintf('\n%-8s %6s %6s | curated: %4s %4s %4s\n', 'mouse', 'cand', ...
    'pool', 'n', 'pos', 'neg');
for iB = 1 : numel(bps)
    [ed, aux] = ed_detect(bps{iB}, 'met', met);
    [~, bname] = fileparts(bps{iB});
    maps = evt_maps(struct('lfp', aux.edSig.lfp), ed.peakTime, aux.fs, ...
        'mapDur', MAPDUR);

    tCur = cur(strcmp({cur.basename}, bname)).peakTime(:);
    isTP = false(numel(ed.peakTime), 1);
    for iC = 1 : numel(tCur)
        [dt, iNear] = min(abs(ed.peakTime - tCur(iC)));
        if dt < TOL, isTP(iNear) = true; end
    end

    keep = evt_gate(ed, met.qa);
    D(iB).bname = bname;
    D(iB).wv    = double(maps.lfp(keep, :));
    D(iB).tst   = maps.tstamps;
    D(iB).tp    = isTP(keep);
    D(iB).nCur  = numel(tCur);
    D(iB).sc    = [ed.fastZ(keep), ed.isoZ(keep), ed.posZ(keep), ...
        ed.amp(keep), ed.dur(keep)];

    % polarity of the CURATED discharges of this mouse
    aCur = ed.amp(isTP);
    fprintf('%-8s %6d %6d | %14d %4d %4d\n', bname(1 : 6), ...
        numel(ed.peakTime), nnz(keep), numel(aCur), nnz(aCur > 0), ...
        nnz(aCur < 0));
end

%% ========================================================================
%  HOW FAR DOES THE RAW PEAK SIT FROM THE FILTERED ONE?
%  ========================================================================
fprintf('\n%-8s %-22s %7s %7s %7s\n', 'mouse', 'set', 'med|dt|', 'dt pos', ...
    'dt neg');
for iB = 1 : numel(D)
    [dt, pol] = peakShift(D(iB).wv, D(iB).tst, SRCH);
    rows = {'pool', true(size(dt)); 'curated', D(iB).tp};
    for iR = 1 : size(rows, 1)
        m = rows{iR, 2};
        fprintf('%-8s %-22s %7.1f %7.1f %7.1f\n', D(iB).bname(1 : 6), ...
            rows{iR, 1}, median(abs(dt(m))) * 1000, ...
            median(dt(m & pol > 0)) * 1000, median(dt(m & pol < 0)) * 1000);
    end
end

%% ========================================================================
%  SWEEP
%  ========================================================================
fprintf('\n%-8s %-6s %-6s %-5s %3s | %5s %4s %4s %6s\n', 'mouse', 'algn', ...
    'detr', 'norm', 'k', 'bestN', 'rec', 'pur', 'load');
res = [];
for iB = 1 : numel(D)
    for iA = 1 : numel(ALGN)
        if strcmp(ALGN{iA}, 'raw')
            [wv, tst] = realign(D(iB).wv, D(iB).tst, SRCH, met.clust.win);
        else
            wv = D(iB).wv; tst = D(iB).tst;
        end
        for iD = 1 : numel(DT)
            for iN = 1 : numel(NM)
                for iK = 1 : numel(KS)
                    cid = ed_clust(wv, tst, 'win', met.clust.win, ...
                        'nPC', met.clust.nPC, 'nClust', KS(iK), ...
                        'scalar', D(iB).sc, 'detrend', DT{iD}, ...
                        'norm', NM{iN});
                    s = scoreClust(cid, D(iB).tp, D(iB).nCur, KS(iK));
                    fprintf(['%-8s %-6s %-6s %-5s %3d | %5d %4.0f %4.0f ' ...
                        '%6d\n'], D(iB).bname(1 : 6), ALGN{iA}, DT{iD}, ...
                        NM{iN}, KS(iK), s.bestN, s.rec * 100, ...
                        s.pur * 100, s.load);
                    res = [res; struct('bname', D(iB).bname, ...
                        'algn', ALGN{iA}, 'detr', DT{iD}, 'norm', NM{iN}, ...
                        'k', KS(iK), 's', s)]; %#ok<AGROW>
                end
            end
        end
    end
end

save(fullfile(DEVDIR, 'ed_alignSweep.mat'), 'res');
fprintf('\nsaved ed_alignSweep.mat\n');

%% ========================================================================
%  SUMMARY: best load over k, per configuration, averaged over mice
%  ========================================================================
fprintf('\n%-6s %-6s %-5s | %s\n', 'algn', 'detr', 'norm', 'mean best load');
cfg = unique(arrayfun(@(r) sprintf('%s|%s|%s', r.algn, r.detr, r.norm), ...
    res, 'uni', false), 'stable');
sm = nan(numel(cfg), 1);
for iC = 1 : numel(cfg)
    p = strsplit(cfg{iC}, '|');
    m = strcmp({res.algn}, p{1}) & strcmp({res.detr}, p{2}) & ...
        strcmp({res.norm}, p{3});
    r = res(m);
    best = nan(1, numel(D));
    for iB = 1 : numel(D)
        l = [r(strcmp({r.bname}, D(iB).bname)).s];
        best(iB) = min([l.load]);
    end
    sm(iC) = mean(best, 'omitnan');
    fprintf('%-6s %-6s %-5s | %8.0f\n', p{1}, p{2}, p{3}, sm(iC));
end
[~, iBest] = min(sm);
fprintf('\nbest: %s\n', cfg{iBest});


% =========================================================================
%  LOCAL
% =========================================================================
function [dt, pol] = peakShift(wv, tst, srch)
% Offset from the detected peak (t = 0) to the raw extremum near it, after the
% slow baseline is taken out - otherwise a wander would set the extremum, not
% the event. POL is the sign of that extremum.
iS = abs(tst) <= srch;
X = wv - mean(wv(:, abs(tst) >= 0.5 * max(abs(tst))), 2);
[~, iPk] = max(abs(X(:, iS)), [], 2);
tS = tst(iS);
dt = tS(iPk)';
pol = sign(arrayfun(@(iEv) X(iEv, find(iS, 1) + iPk(iEv) - 1), ...
    (1 : size(X, 1))'));

end     % peakShift


function [wv, tst] = realign(wv, tst, srch, win)
% Re-cut each waveform so t = 0 sits at its own raw extremum. The map is wider
% than the clustering window, so this is an index shift, not a re-read.
dt = peakShift(wv, tst, srch);
fs = 1 / median(diff(tst));
sh = round(dt * fs);
keep = abs(tst) <= max(abs(win)) + srch;
i0 = find(keep);
out = nan(size(wv, 1), numel(i0));
for iEv = 1 : size(wv, 1)
    idx = i0 + sh(iEv);
    ok = idx >= 1 & idx <= size(wv, 2);
    out(iEv, ok) = wv(iEv, idx(ok));
end
wv = out;
tst = tst(i0) - median(tst(i0));

end     % realign


function s = scoreClust(cid, tp, nCur, k)
s = struct('bestN', 0, 'rec', 0, 'pur', 0, 'load', NaN);
dens = zeros(k, 1); cnt = zeros(k, 1); hit = zeros(k, 1);
for iK = 1 : k
    inK = cid == iK;
    cnt(iK) = nnz(inK);
    hit(iK) = nnz(inK & tp);
    if cnt(iK) > 0, dens(iK) = hit(iK) / cnt(iK); end
    r = hit(iK) / max(1, nCur);
    if r > s.rec
        s.rec = r; s.bestN = cnt(iK);
        s.pur = hit(iK) / max(1, cnt(iK));
    end
end
[~, ord] = sort(dens, 'descend');
cumHit = cumsum(hit(ord)); cumN = cumsum(cnt(ord));
iEnough = find(cumHit >= 0.8 * nCur, 1);
if ~isempty(iEnough), s.load = cumN(iEnough); end

end     % scoreClust
