% ED_CLUSTSWEEP  Can waveform clustering isolate the discharges?
%
% The question that decides whether ed_curate can replace per-event curation.
% Sweeps what the clustering is given:
%   win  - waveform window. Short isolates the spike; long also sees the decay,
%          where a discharge (back to baseline in 50-100 ms) and a sharp wave
%          (300+ ms) differ most.
%   norm - raw amplitude, or each waveform scaled to unit peak. Amplitude
%          separates discharges from ordinary deflections, but it also lets the
%          few largest events dominate the principal components.
%   feat - waveform components alone, or plus the scalar shape measures already
%          computed per event (fastZ isoZ posZ amp dur).
%   k    - cluster count.
%
% Scored, for the cluster richest in curated discharges:
%   rec  - fraction of curated discharges it holds
%   pur  - fraction of it that is curated discharges
%   load - events to review to reach 80% recall, taking clusters in order of
%          discharge density. THE number the user feels; compare it to the
%          pool size, which is what per-event curation would cost.
%
% Writes ed_clustSweep.mat.

clear
DEVDIR = fileparts(mfilename('fullpath'));
TOL = 0.050;
MAPDUR = [-0.1 0.1];

bps = { ...
    'D:\Data\RA\raMCU3\raMCU3_211203_084720', ...
    'D:\Data\RA\raMCU4\raMCU4_211220_0834', ...
    'D:\Data\RA\raMCU5\raMCU5_220322_1906'};

WINS  = {[-0.015 0.015], [-0.025 0.025], [-0.05 0.05], [-0.1 0.1]};
NORMS = {'raw', 'unit'};
FEATS = {'pca', 'pca+scalar'};
KS    = [4 6 8 12];

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
    D(iB).bname  = bname;
    D(iB).wv     = double(maps.lfp(keep, :));
    D(iB).tst    = maps.tstamps;
    D(iB).tp     = isTP(keep);
    D(iB).nCur   = numel(tCur);
    D(iB).scalar = [ed.fastZ(keep), ed.isoZ(keep), ed.posZ(keep), ...
        ed.amp(keep), ed.dur(keep)];
    fprintf('%s: ch%d, %d cand, %d pooled, %d/%d curated inside\n', ...
        bname, ed.info.edCh, numel(ed.peakTime), nnz(keep), ...
        nnz(isTP & keep), numel(tCur));
end

fprintf('\n%-7s %-11s %-5s %-11s %3s | %5s %4s %4s %6s %6s\n', 'mouse', ...
    'win(ms)', 'norm', 'feat', 'k', 'bestN', 'rec', 'pur', 'load', 'pool');
res = [];
for iB = 1 : numel(D)
    for iW = 1 : numel(WINS)
        for iN = 1 : numel(NORMS)
            W = D(iB).wv;
            if strcmp(NORMS{iN}, 'unit')
                W = W ./ max(abs(W), [], 2);
            end
            for iF = 1 : numel(FEATS)
                if strcmp(FEATS{iF}, 'pca'), sc = []; else
                    sc = D(iB).scalar;
                end
                for iK = 1 : numel(KS)
                    cid = ed_clust(W, D(iB).tst, 'win', WINS{iW}, ...
                        'nPC', 6, 'nClust', KS(iK), 'scalar', sc);
                    s = scoreClust(cid, D(iB).tp, D(iB).nCur, KS(iK));
                    fprintf(['%-7s %-11s %-5s %-11s %3d | %5d %4.0f %4.0f ' ...
                        '%6d %6d\n'], D(iB).bname(1 : 6), ...
                        mat2str(WINS{iW} * 1000), NORMS{iN}, FEATS{iF}, ...
                        KS(iK), s.bestN, s.rec * 100, s.pur * 100, ...
                        s.load, numel(D(iB).tp));
                    res = [res; struct('bname', D(iB).bname, ...
                        'win', WINS{iW}, 'norm', NORMS{iN}, ...
                        'feat', FEATS{iF}, 'k', KS(iK), 's', s, ...
                        'pool', numel(D(iB).tp))]; %#ok<AGROW>
                end
            end
        end
    end
end

save(fullfile(DEVDIR, 'ed_clustSweep.mat'), 'res');
fprintf('\nsaved ed_clustSweep.mat\n');


% =========================================================================
%  LOCAL
% =========================================================================
function s = scoreClust(cid, tp, nCur, k)
% Best cluster + the review load to reach 80% recall over ranked clusters.
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
