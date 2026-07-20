% ED_SIGSWEEP  Which LFP signal should ED detection run on?
%
% Moving from the sleep_sig average to one auto-picked channel made detection
% much noisier (raMCU5: 5131 candidates before, 13457 after; 45 gated before,
% 1254 after). That is a real finding and it points somewhere I dismissed too
% fast: a discharge appears on EVERY channel at once (measured spread 0.6-0.9
% in ed_laminar.m), so averaging channels is a matched filter for it - it
% suppresses uncorrelated noise and keeps the event. The fault in sleep_sig was
% not that it averaged, but that it averaged an arbitrary, per-mouse set.
%
% So compare, on the same detector and the same curated ground truth:
%   pick   - the single ed_pickCh channel
%   best   - the single channel where the curated discharges are largest
%   allCh  - the mean of every neural channel (fully deterministic, no pick)
%   grp    - the mean of the spike group holding the pick
%
% The number that decides it is AUC: how well fastZ separates the curated
% discharges from every other candidate. Recall is not enough - a signal can
% propose everything and separate nothing.
%
% Writes ed_sigSweep.mat.

clear
DEVDIR = fileparts(mfilename('fullpath'));
TOL = 0.050;

bps = { ...
    'D:\Data\RA\raMCU3\raMCU3_211203_084720', ...
    'D:\Data\RA\raMCU4\raMCU4_211220_0834', ...
    'D:\Data\RA\raMCU5\raMCU5_220322_1906'};
bestCh = [13, 7, 9];        % per-channel edAmp winner, from ed_laminar.m

cur = load(fullfile(DEVDIR, 'ed_curatedTimes.mat'));
cur = cur.edCurated;
met = ed_methods('default');

res = struct([]);
for iB = 1 : numel(bps)

    basepath = bps{iB};
    [~, basename] = fileparts(basepath);
    v = basepaths2vars('basepaths', {basepath}, 'vars', {'session'}, ...
        'flgPrnt', false);
    ex = v.session.extracellular;
    chAll = sort(unique([ex.spikeGroups.channels{:}]));

    chPick = ed_pickCh(basepath, 'flgForce', true);
    iGrp = find(cellfun(@(c) ismember(chPick, c), ...
        ex.spikeGroups.channels), 1);

    opt = struct( ...
        'pick',  chPick, ...
        'best',  bestCh(iB), ...
        'allCh', chAll, ...
        'grp',   ex.spikeGroups.channels{iGrp}(:)');

    tCur = cur(strcmp({cur.basename}, basename)).peakTime(:);
    fprintf('\n===== %s (%d curated) =====\n', basename, numel(tCur));
    fprintf('%-7s %7s %7s %7s %7s %7s\n', 'signal', 'ch', 'nCand', ...
        'recall', 'AUC', 'nGate');

    fn = fieldnames(opt);
    for iO = 1 : numel(fn)
        [ed, ~] = ed_detect(basepath, 'met', met, 'edCh', opt.(fn{iO}));

        isTP = false(numel(ed.peakTime), 1);
        for iC = 1 : numel(tCur)
            [dt, iNear] = min(abs(ed.peakTime - tCur(iC)));
            if dt < TOL, isTP(iNear) = true; end
        end
        rec = nnz(isTP) / numel(tCur);
        auc = rocAuc(ed.fastZ(isTP), ed.fastZ(~isTP));
        nGate = nnz(evt_gate(ed, met.qa));

        fprintf('%-7s %7s %7d %6.0f%% %7.3f %7d\n', fn{iO}, ...
            mat2str(opt.(fn{iO})(1 : min(2, end))), numel(ed.peakTime), ...
            rec * 100, auc, nGate);

        res(iB).(fn{iO}) = struct('ch', opt.(fn{iO}), ...
            'nCand', numel(ed.peakTime), 'rec', rec, 'auc', auc, ...
            'nGate', nGate);
    end
    res(iB).basename = basename;
end

save(fullfile(DEVDIR, 'ed_sigSweep.mat'), 'res');
fprintf('\nsaved ed_sigSweep.mat\n');


% =========================================================================
%  LOCAL
% =========================================================================
function a = rocAuc(pos, neg)
% Area under the ROC via the rank-sum identity.
pos = pos(isfinite(pos)); neg = neg(isfinite(neg));
if isempty(pos) || isempty(neg), a = NaN; return, end
r = tiedrank([pos(:); neg(:)]);
a = (sum(r(1 : numel(pos))) - numel(pos) * (numel(pos) + 1) / 2) ...
    / (numel(pos) * numel(neg));
end
