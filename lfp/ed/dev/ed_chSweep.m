% ED_CHSWEEP  Per-channel ground truth for the channel picker.
%
% Runs the real detector on EVERY neural channel of raMCU3/4/5 and scores each
% by the number that matters: the AUC with which fastZ separates the curated
% discharges from every other candidate on that channel. Then asks which cheap,
% detection-free criterion predicts it:
%
%   p99.9  - prctile(|filt|, 99.9)          (what ed_pickCh does now)
%   xRate  - fraction of samples with |filt| > thr * robust scale
%
% xRate is the principled candidate: every channel sees the discharge (measured
% spread 0.6-0.9 in ed_laminar.m), so the best channel is simply the one that
% produces the FEWEST false candidates - the least noise in the detection band.
%
% Writes ed_chSweep.mat. Slow (one full detection per channel).

clear
DEVDIR = fileparts(mfilename('fullpath'));
TOL = 0.050;
NPROBE = 20; PROBDUR = 15;

bps = { ...
    'D:\Data\RA\raMCU3\raMCU3_211203_084720', ...
    'D:\Data\RA\raMCU4\raMCU4_211220_0834', ...
    'D:\Data\RA\raMCU5\raMCU5_220322_1906'};

cur = load(fullfile(DEVDIR, 'ed_curatedTimes.mat'));
cur = cur.edCurated;
met = ed_methods('default');

res = cell(numel(bps), 1);
for iB = 1 : numel(bps)

    basepath = bps{iB};
    [~, basename] = fileparts(basepath);
    v = basepaths2vars('basepaths', {basepath}, 'vars', {'session'}, ...
        'flgPrnt', false);
    ex = v.session.extracellular;
    fs = ex.srLfp; nCh = ex.nChannels;
    chOk = sort(unique([ex.spikeGroups.channels{:}]));
    if round(ex.sr) == 24414, b2u = 1; else, b2u = 0.195; end
    tCur = cur(strcmp({cur.basename}, basename)).peakTime(:);

    % ------------------------------------------- cheap criteria, probe data
    fname = fullfile(basepath, [basename '.lfp']);
    d = dir(fname);
    durRec = d.bytes / 2 / nCh / fs;
    tProbe = linspace(0, durRec - PROBDUR, NPROBE);
    sig = cell(NPROBE, 1);
    for iPrb = 1 : NPROBE
        sig{iPrb} = double(binary_load(fname, 'fs', fs, 'nCh', nCh, ...
            'start', tProbe(iPrb), 'duration', PROBDUR, 'ch', chOk, ...
            'bit2uv', b2u));
    end
    sig = vertcat(sig{:});

    pk = zeros(1, numel(chOk)); xRate = pk;
    for iCh = 1 : numel(chOk)
        fb = filterLFP(sig(:, iCh), 'fs', fs, 'type', 'butter', ...
            'dataOnly', true, 'order', 3, 'passband', met.passband, ...
            'graphics', false);
        scl = 1.4826 * median(abs(fb - median(fb)));
        pk(iCh) = prctile(abs(fb), 99.9);
        xRate(iCh) = mean(abs(fb) > met.thr * scl);
    end

    % -------------------------------------------------- ground truth per ch
    auc = nan(1, numel(chOk)); nCand = auc; nGate = auc;
    for iCh = 1 : numel(chOk)
        ed = ed_detect(basepath, 'met', met, 'edCh', chOk(iCh));
        isTP = false(numel(ed.peakTime), 1);
        for iC = 1 : numel(tCur)
            [dt, iNear] = min(abs(ed.peakTime - tCur(iC)));
            if dt < TOL, isTP(iNear) = true; end
        end
        auc(iCh) = rocAuc(ed.fastZ(isTP), ed.fastZ(~isTP));
        nCand(iCh) = numel(ed.peakTime);
        nGate(iCh) = nnz(evt_gate(ed, met.qa));
    end

    fprintf('\n===== %s =====\n', basename);
    fprintf('%4s %10s %10s %8s %8s %8s\n', 'ch', 'p99.9', 'xRate', ...
        'nCand', 'nGate', 'AUC');
    for iCh = 1 : numel(chOk)
        fprintf('%4d %10.1f %10.5f %8d %8d %8.3f\n', chOk(iCh), pk(iCh), ...
            xRate(iCh), nCand(iCh), nGate(iCh), auc(iCh));
    end
    [~, iP] = max(pk); [~, iX] = min(xRate); [~, iA] = max(auc);
    fprintf('  p99.9 picks ch%d (AUC %.3f) | xRate picks ch%d (AUC %.3f)', ...
        chOk(iP), auc(iP), chOk(iX), auc(iX));
    fprintf(' | best is ch%d (AUC %.3f)\n', chOk(iA), auc(iA));
    fprintf('  corr(p99.9,AUC) = %.2f | corr(-xRate,AUC) = %.2f\n', ...
        corr(pk(:), auc(:)), corr(-xRate(:), auc(:)));

    res{iB} = struct('basename', basename, 'ch', chOk, 'pk', pk, ...
        'xRate', xRate, 'nCand', nCand, 'nGate', nGate, 'auc', auc);
end

save(fullfile(DEVDIR, 'ed_chSweep.mat'), 'res');
fprintf('\nsaved ed_chSweep.mat\n');


% =========================================================================
%  LOCAL
% =========================================================================
function a = rocAuc(pos, neg)
pos = pos(isfinite(pos)); neg = neg(isfinite(neg));
if isempty(pos) || isempty(neg), a = NaN; return, end
r = tiedrank([pos(:); neg(:)]);
a = (sum(r(1 : numel(pos))) - numel(pos) * (numel(pos) + 1) / 2) ...
    / (numel(pos) * numel(neg));
end
