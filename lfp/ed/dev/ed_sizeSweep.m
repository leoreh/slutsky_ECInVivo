% ED_SIZESWEEP  How much should SIZE count in the clustering?
%
% Leore's objection: amplitude carries most of the variance, he uses it when
% curating by hand, and the clustering puts a 3400 uV event with 1000 uV ones.
% Unit-peak normalisation removes size from the waveform on purpose - otherwise
% the components describe nothing else - and the scalar block hands back only
% its RANK, which keeps the order and throws away the ratio. So no setting of
% the current parameters can express "three times larger".
%
% evt_clust now takes wSize: log10(size) as one explicit axis, scaled to the
% spread of the leading shape component. wSize = 0 is shape-only (the shipped
% behaviour), 1 means size counts as much as the dominant shape axis.
%
% Two things must be traded off, so both are measured:
%   HOMOGENEITY - do clusters stop mixing amplitudes? Median over clusters of
%     the p90/p10 amplitude ratio within a cluster. This is Leore's complaint,
%     measured directly.
%   RECALL - does the discharge still concentrate? Best-cluster recall of the
%     curated events, and the review load to reach 80% of them.

clear
DEVDIR = fileparts(mfilename('fullpath'));
TOL = 0.050;
MAPDUR = [-0.1 0.1];
WS = [0 0.5 1 2 4];

bps = { ...
    'D:\Data\RA\raMCU3\raMCU3_211203_084720', ...
    'D:\Data\RA\raMCU4\raMCU4_211220_0834', ...
    'D:\Data\RA\raMCU5\raMCU5_220322_1906'};

met = ed_methods('default');
cur = load(fullfile(DEVDIR, 'ed_curatedTimes.mat'));
cur = cur.edCurated;

D = struct([]);
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
    D(iB).wv   = double(maps.lfp(keep, :));
    D(iB).tst  = maps.tstamps;
    D(iB).tp   = isTP(keep);
    D(iB).nCur = numel(tCur);
    D(iB).amp  = abs(ed.amp(keep));
    D(iB).sc   = [ed.fastZ(keep), ed.isoZ(keep), ed.posZ(keep), ...
        ed.amp(keep), ed.dur(keep)];
end

fprintf('\n%-8s %6s | %6s %6s %6s %6s\n', 'mouse', 'wSize', 'ampRat', ...
    'rec', 'pur', 'load');
res = [];
for iB = 1 : numel(D)
    for iW = 1 : numel(WS)
        cid = evt_clust(D(iB).wv, D(iB).tst, 'win', met.clust.win, ...
            'nPC', met.clust.nPC, 'nClust', met.clust.nClust, ...
            'scalar', D(iB).sc, 'detrend', met.clust.detrend, ...
            'norm', met.clust.norm, 'wSize', WS(iW));
        s = scoreClust(cid, D(iB).tp, D(iB).nCur);
        s.ampRat = ampSpread(cid, D(iB).amp);
        fprintf('%-8s %6.1f | %6.1f %6.0f %6.0f %6d\n', D(iB).bname(1 : 6), ...
            WS(iW), s.ampRat, s.rec * 100, s.pur * 100, s.load);
        res = [res; struct('bname', D(iB).bname, 'w', WS(iW), 's', s)]; %#ok
    end
end

fprintf('\n%-6s | %8s %8s %8s\n', 'wSize', 'ampRat', 'rec %', 'load');
for iW = 1 : numel(WS)
    m = [res.w] == WS(iW);
    r = [res(m).s];
    fprintf('%-6.1f | %8.1f %8.0f %8.0f\n', WS(iW), ...
        mean([r.ampRat]), mean([r.rec]) * 100, mean([r.load]));
end

save(fullfile(DEVDIR, 'ed_sizeSweep.mat'), 'res');
fprintf('\nsaved ed_sizeSweep.mat\n');


% =========================================================================
%  LOCAL
% =========================================================================
function r = ampSpread(cid, amp)
% Median over clusters of the within-cluster p90/p10 amplitude ratio. A cluster
% holding a 3400 and a 1000 scores 3.4; one holding only its own kind scores
% near 1.
k = max(cid);
v = nan(k, 1);
for iK = 1 : k
    a = amp(cid == iK);
    if numel(a) < 5, continue; end
    q = prctile(a, [10 90]);
    if q(1) > 0, v(iK) = q(2) / q(1); end
end
r = median(v, 'omitnan');

end     % ampSpread


function s = scoreClust(cid, tp, nCur)
s = struct('bestN', 0, 'rec', 0, 'pur', 0, 'load', NaN);
k = max(cid);
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
