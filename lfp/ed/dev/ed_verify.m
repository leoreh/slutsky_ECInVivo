% ED_VERIFY  End-to-end check of the rebuilt pipeline across the cohort.
%
% Runs detect -> noise filter -> cluster on every MCU session, in memory,
% writing nothing. Reports per session: the auto-picked channel, candidates,
% pool size, clusters, and - where hand-curated discharges exist - recall and
% the review load the clustering implies.
%
% What to look for: the CAG mice (raMCU*) should carry discharges and the
% controls (lh*) should not, so the pool on a control is junk the user rejects
% wholesale. Recall is the number that must not regress.

clear
DEVDIR = fileparts(mfilename('fullpath'));
TOL = 0.050;
MAPDUR = [-0.1 0.1];

bps = [mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl'), ...
    mcu_basepaths('ra')];

cur = load(fullfile(DEVDIR, 'ed_curatedTimes.mat'));
cur = cur.edCurated;
met = ed_methods('default');

fprintf('%-22s %4s %6s %7s %6s %5s %6s %6s\n', 'session', 'ch', 'h', ...
    'nCand', 'pool', 'nCl', 'recall', 'review');
rows = cell(numel(bps), 1);
for iB = 1 : numel(bps)

    basepath = bps{iB};
    [~, basename] = fileparts(basepath);
    t0 = tic;

    [ed, aux] = ed_detect(basepath, 'met', met);
    pool = evt_gate(ed, met.qa);
    maps = evt_maps(struct('lfp', aux.edSig.lfp), ed.peakTime, aux.fs, ...
        'mapDur', MAPDUR);

    iPool = find(pool);
    cid = ed_clust(maps.lfp(iPool, :), maps.tstamps, 'win', met.clust.win, ...
        'nPC', met.clust.nPC, 'nClust', met.clust.nClust, 'scalar', ...
        [ed.fastZ(iPool), ed.isoZ(iPool), ed.posZ(iPool), ...
        ed.amp(iPool), ed.dur(iPool)]);

    rec = NaN; load80 = NaN;
    iCur = find(strcmp({cur.basename}, basename));
    if ~isempty(iCur)
        tCur = cur(iCur).peakTime(:);
        isTP = false(numel(ed.peakTime), 1);
        for iC = 1 : numel(tCur)
            [dt, iNear] = min(abs(ed.peakTime - tCur(iC)));
            if dt < TOL, isTP(iNear) = true; end
        end
        rec = nnz(isTP & pool) / numel(tCur);
        load80 = reviewLoad(cid, isTP(iPool), numel(tCur));
    end

    fprintf('%-22s %4d %6.1f %7d %6d %5d %5.0f%% %6d   (%.0fs)\n', ...
        basename, ed.info.edCh, ed.info.sigDur / 3600, ...
        numel(ed.peakTime), nnz(pool), max(cid), rec * 100, load80, toc(t0));

    rows{iB} = struct('basename', basename, 'edCh', ed.info.edCh, ...
        'nCand', numel(ed.peakTime), 'pool', nnz(pool), 'rec', rec, ...
        'load80', load80);
end

res = [rows{:}];
save(fullfile(DEVDIR, 'ed_verify.mat'), 'res');
fprintf('\nsaved ed_verify.mat\n');


% =========================================================================
%  LOCAL
% =========================================================================
function n = reviewLoad(cid, tp, nCur)
% Events to review to reach 80% of the curated discharges, taking clusters in
% order of discharge density.
n = NaN;
k = max(cid);
if isempty(k) || isnan(k), return, end
cnt = zeros(k, 1); hit = zeros(k, 1);
for iK = 1 : k
    inK = cid == iK;
    cnt(iK) = nnz(inK);
    hit(iK) = nnz(inK & tp);
end
dens = hit ./ max(1, cnt);
[~, ord] = sort(dens, 'descend');
iEnough = find(cumsum(hit(ord)) >= 0.8 * nCur, 1);
if ~isempty(iEnough)
    cumN = cumsum(cnt(ord));
    n = cumN(iEnough);
end
end
