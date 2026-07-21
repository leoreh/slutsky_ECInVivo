% ED_GATESWEEP  How should the noise filter hand events to the clustering?
%
% On the single raw channel a FIXED fastZ/isoZ threshold does not transfer:
% the same numbers leave 292 events in raMCU3 and 1254 in raMCU5, because the
% two recordings differ in background, not in discharge count. So compare:
%   (a) fixed z thresholds, swept
%   (b) a RANK gate - the top N candidates by fastZ, which fixes the curation
%       load by construction
% and for each, measure recall of the curated discharges AND whether the
% clustering can then isolate them.
%
% Writes ed_gateSweep.mat.

clear
DEVDIR = fileparts(mfilename('fullpath'));
TOL = 0.050;
MAPDUR = [-0.1 0.1];
FASTZ = [15 20 25 30 40 50];
TOPN  = [100 200 300 500];

bps = { ...
    'D:\Data\RA\raMCU3\raMCU3_211203_084720', ...
    'D:\Data\RA\raMCU4\raMCU4_211220_0834', ...
    'D:\Data\RA\raMCU5\raMCU5_220322_1906'};

cur = load(fullfile(DEVDIR, 'ed_curatedTimes.mat'));
cur = cur.edCurated;
met = ed_methods('default');

res = struct([]);
for iB = 1 : numel(bps)

    basepath = bps{iB};
    [~, basename] = fileparts(basepath);

    [ed, aux] = ed_detect(basepath, 'met', met);
    maps = evt_maps(struct('lfp', aux.edSig.lfp), ed.peakTime, aux.fs, ...
        'mapDur', MAPDUR);

    tCur = cur(strcmp({cur.basename}, basename)).peakTime(:);
    isTP = false(numel(ed.peakTime), 1);
    for iC = 1 : numel(tCur)
        [dt, iNear] = min(abs(ed.peakTime - tCur(iC)));
        if dt < TOL, isTP(iNear) = true; end
    end

    fprintf('\n===== %s : %d cand, %d curated =====\n', basename, ...
        numel(ed.peakTime), numel(tCur));

    % --------------------------------------------------- (a) fixed threshold
    fprintf('  fixed fastZ (isoZ >= 20):\n');
    for iF = 1 : numel(FASTZ)
        k = ed.fastZ >= FASTZ(iF) & ed.isoZ >= 20;
        fprintf('    fastZ>=%2d : %5d kept, recall %3.0f%%\n', FASTZ(iF), ...
            nnz(k), nnz(k & isTP) / numel(tCur) * 100);
    end

    % ------------------------------------------------------- (b) rank gate
    fprintf('  rank gate (top N by fastZ, isoZ >= 20 first):\n');
    for iN = 1 : numel(TOPN)
        ok = find(ed.isoZ >= 20);
        [~, ord] = sort(ed.fastZ(ok), 'descend');
        iKeep = ok(ord(1 : min(TOPN(iN), numel(ord))));
        rec = nnz(isTP(iKeep)) / numel(tCur);

        [cid, cInfo] = evt_clust(maps.lfp(iKeep, :), maps.tstamps, ...
            'win', met.clust.win, 'nPC', met.clust.nPC, ...
            'nClust', met.clust.nClust);
        tpIn = isTP(iKeep);

        bRec = 0; bPur = 0; bN = 0;
        for iK = 1 : cInfo.nClust
            inK = cid == iK;
            r = nnz(inK & tpIn) / max(1, numel(tCur));
            if r > bRec
                bRec = r; bN = nnz(inK);
                bPur = nnz(inK & tpIn) / max(1, nnz(inK));
            end
        end
        nSpread = numel(unique(cid(tpIn & ~isnan(cid))));

        fprintf(['    N=%3d : recall %3.0f%%, %d clusters | best n=%3d ' ...
            'rec %3.0f%% pur %3.0f%% | spread %d\n'], TOPN(iN), rec * 100, ...
            cInfo.nClust, bN, bRec * 100, bPur * 100, nSpread);

        res(iB).topN(iN) = struct('n', TOPN(iN), 'rec', rec, ...
            'nClust', cInfo.nClust, 'bestN', bN, 'bestRec', bRec, ...
            'bestPur', bPur, 'spread', nSpread);
    end

    res(iB).basename = basename;
    res(iB).isTP = isTP;
    res(iB).fastZ = ed.fastZ;
    res(iB).isoZ = ed.isoZ;
end

save(fullfile(DEVDIR, 'ed_gateSweep.mat'), 'res');
fprintf('\nsaved ed_gateSweep.mat\n');
