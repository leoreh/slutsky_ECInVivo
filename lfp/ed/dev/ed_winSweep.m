% ED_WINSWEEP  Is the clustering window too wide, at the CURRENT filter?
%
% ed_clustSweep.m chose +-50 ms, but it ran under the old gate (fastZ >= 15,
% isoZ >= 20) where a pool was 75-490 events. The shipped filter is now much
% more permissive and a pool is 2700-8500, so the earlier answer does not
% automatically hold: with a wider pool the window carries proportionally more
% baseline, and baseline is what the principal components would then describe.
%
% Also tests whether AMPLITUDE is doing the work. The clusters in a real
% session look like one shape at a dozen sizes, which is what you would see if
% the scalar measures (which include .amp) dominated the geometry rather than
% the waveform.
%
% Scored, as before, against the curated discharges: events you must review to
% reach 80% of them, taking clusters in order of discharge density.

clear
DEVDIR = fileparts(mfilename('fullpath'));
TOL = 0.050;
MAPDUR = [-0.1 0.1];

bps = { ...
    'D:\Data\RA\raMCU3\raMCU3_211203_084720', ...
    'D:\Data\RA\raMCU4\raMCU4_211220_0834', ...
    'D:\Data\RA\raMCU5\raMCU5_220322_1906'};

WINS = {[-0.05 0.05]};
SCAL = {'shape', 'all'};           % which scalar measures join the PCs
KS   = [12 20 30 45 60];

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
    D(iB).wv    = double(maps.lfp(keep, :));
    D(iB).tst   = maps.tstamps;
    D(iB).tp    = isTP(keep);
    D(iB).nCur  = numel(tCur);
    % 'shape' drops .amp, which is the one measure that is purely size
    D(iB).sAll   = [ed.fastZ(keep), ed.isoZ(keep), ed.posZ(keep), ...
        ed.amp(keep), ed.dur(keep)];
    D(iB).sShape = [ed.fastZ(keep), ed.isoZ(keep), ed.dur(keep)];
    fprintf('%s: %d cand, %d pooled, %d/%d curated inside\n', bname, ...
        numel(ed.peakTime), nnz(keep), nnz(isTP & keep), numel(tCur));
end

fprintf('\n%-7s %-9s %-6s %3s | %5s %4s %4s %6s %6s\n', 'mouse', 'win(ms)', ...
    'scalar', 'k', 'bestN', 'rec', 'pur', 'load', 'pool');
res = [];
for iB = 1 : numel(D)
    for iW = 1 : numel(WINS)
        for iS = 1 : numel(SCAL)
            switch SCAL{iS}
                case 'none',  sc = [];
                case 'shape', sc = D(iB).sShape;
                otherwise,    sc = D(iB).sAll;
            end
            for iK = 1 : numel(KS)
                cid = ed_clust(D(iB).wv, D(iB).tst, 'win', WINS{iW}, ...
                    'nPC', met.clust.nPC, 'nClust', KS(iK), 'scalar', sc);
                s = scoreClust(cid, D(iB).tp, D(iB).nCur, KS(iK));
                fprintf('%-7s %-9s %-6s %3d | %5d %4.0f %4.0f %6d %6d\n', ...
                    D(iB).bname(1 : 6), mat2str(WINS{iW} * 1000), ...
                    SCAL{iS}, KS(iK), s.bestN, s.rec * 100, s.pur * 100, ...
                    s.load, numel(D(iB).tp));
                res = [res; struct('bname', D(iB).bname, 'win', WINS{iW}, ...
                    'scalar', SCAL{iS}, 'k', KS(iK), 's', s, ...
                    'pool', numel(D(iB).tp))]; %#ok<AGROW>
            end
        end
    end
end

save(fullfile(DEVDIR, 'ed_winSweep.mat'), 'res');
fprintf('\nsaved ed_winSweep.mat\n');


% =========================================================================
%  LOCAL
% =========================================================================
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
end
