% ED_CLUSTTEST  The make-or-break test for cluster curation.
%
% Two questions, both against the 47 hand-curated discharges in raMCU3/4/5:
%   1. RECALL. Detection now runs on one auto-picked .lfp channel instead of
%      the sleep_sig average. Are the curated discharges still proposed, and
%      do they still survive the noise filter?
%   2. SEPARABILITY. Does PCA + GMM put them in ONE cluster, and is that
%      cluster mostly them? If a discharge type does not come out as a cluster,
%      cluster curation cannot replace per-event curation and I have to say so.
%
% Reports, per mouse: candidates, gated survivors, recall at each stage, and
% for the best-matching cluster its size, recall and purity.
%
% Writes ed_clustTest.mat.

clear
DEVDIR = fileparts(mfilename('fullpath'));
TOL = 0.050;                % a curated peak is matched within this [s]

bps = { ...
    'D:\Data\RA\raMCU1\raMCU1_080621_0930', ...
    'D:\Data\RA\raMCU2\raMCU2_080621_0930', ...
    'D:\Data\RA\raMCU3\raMCU3_211203_084720', ...
    'D:\Data\RA\raMCU4\raMCU4_211220_0834', ...
    'D:\Data\RA\raMCU5\raMCU5_220322_1906', ...
    'D:\Data\lh100\lh100_220413_111004', ...
    'D:\Data\lh107\lh107_220518_091200', ...
    'D:\Data\lh132\lh132_230413_094013'};

cur = load(fullfile(DEVDIR, 'ed_curatedTimes.mat'));
cur = cur.edCurated;
met = ed_methods('default');
MAPDUR = [-0.1 0.1];

res = struct([]);
for iB = 1 : numel(bps)

    basepath = bps{iB};
    [~, basename] = fileparts(basepath);

    [ed, aux] = ed_detect(basepath, 'met', met);
    nCand = numel(ed.peakTime);
    keep = evt_gate(ed, met.qa);

    maps = evt_maps(struct('lfp', aux.edSig.lfp), ed.peakTime, aux.fs, ...
        'mapDur', MAPDUR);

    % cluster the survivors only
    iKeep = find(keep);
    [cid, cInfo] = ed_clust(maps.lfp(iKeep, :), maps.tstamps, ...
        'win', met.clust.win, 'nPC', met.clust.nPC, ...
        'nClust', met.clust.nClust, 'scalar', ...
        [ed.fastZ(iKeep), ed.isoZ(iKeep), ed.posZ(iKeep), ...
        ed.amp(iKeep), ed.dur(iKeep)]);

    r = struct('basename', basename, 'edCh', ed.info.edCh, ...
        'nCand', nCand, 'nKeep', nnz(keep), 'nClust', cInfo.nClust, ...
        'durH', ed.info.sigDur / 3600);

    % ---------------------------------------------------------- ground truth
    iCur = find(strcmp({cur.basename}, basename));
    if ~isempty(iCur)
        tCur = cur(iCur).peakTime(:);
        isTP = false(nCand, 1);
        for iC = 1 : numel(tCur)
            [dt, iNear] = min(abs(ed.peakTime - tCur(iC)));
            if dt < TOL, isTP(iNear) = true; end
        end
        r.nCur    = numel(tCur);
        r.recDet  = nnz(isTP) / numel(tCur);
        r.recGate = nnz(isTP & keep) / numel(tCur);

        % best cluster for the curated discharges
        tpIn = isTP(iKeep);
        best = struct('id', NaN, 'n', 0, 'rec', 0, 'pur', 0);
        for iK = 1 : cInfo.nClust
            inK = cid == iK;
            rec = nnz(inK & tpIn) / max(1, numel(tCur));
            if rec > best.rec
                best = struct('id', iK, 'n', nnz(inK), 'rec', rec, ...
                    'pur', nnz(inK & tpIn) / max(1, nnz(inK)));
            end
        end
        r.best = best;
        % how many clusters hold at least one curated discharge
        r.nClustWithTP = numel(unique(cid(tpIn & ~isnan(cid))));
    end

    res(iB).r = r;
    res(iB).cid = cid;
    res(iB).cInfo = cInfo;

    fprintf('\n%s (ch %d, %.1f h)\n', basename, r.edCh, r.durH);
    fprintf('  %d cand -> %d gated -> %d clusters\n', r.nCand, r.nKeep, ...
        r.nClust);
    if isfield(r, 'nCur')
        fprintf('  curated %d : recall det %.0f%%, gated %.0f%%\n', ...
            r.nCur, r.recDet * 100, r.recGate * 100);
        fprintf(['  best cluster #%d: n=%d, holds %.0f%% of curated, ', ...
            '%.0f%% pure | curated spread over %d clusters\n'], ...
            r.best.id, r.best.n, r.best.rec * 100, r.best.pur * 100, ...
            r.nClustWithTP);
    end
end

save(fullfile(DEVDIR, 'ed_clustTest.mat'), 'res');
fprintf('\nsaved ed_clustTest.mat\n');
