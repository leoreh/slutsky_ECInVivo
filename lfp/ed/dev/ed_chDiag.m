% ED_CHDIAG  Per-channel breakdown of the ed_pickCh score against the
% per-channel curated-discharge amplitude, to see WHY a pick goes wrong.

DEVDIR = fileparts(mfilename('fullpath'));
load(fullfile(DEVDIR, 'ed_laminar.mat'), 'res');

bps = { ...
    'D:\Data\RA\raMCU1\raMCU1_080621_0930', ...
    'D:\Data\RA\raMCU2\raMCU2_080621_0930', ...
    'D:\Data\RA\raMCU3\raMCU3_211203_084720', ...
    'D:\Data\RA\raMCU4\raMCU4_211220_0834', ...
    'D:\Data\RA\raMCU5\raMCU5_220322_1906'};

NPROBE = 20; PROBDUR = 15; PASSBAND = [60 150];

for iB = 1 : numel(bps)
    basepath = bps{iB};
    [~, basename] = fileparts(basepath);

    v = basepaths2vars('basepaths', {basepath}, 'vars', {'session'}, ...
        'flgPrnt', false);
    ex = v.session.extracellular;
    fs = ex.srLfp; nCh = ex.nChannels;
    chOk = sort(unique([ex.spikeGroups.channels{:}]));
    if round(ex.sr) == 24414, b2u = 1; else, b2u = 0.195; end

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

    scl = zeros(1, numel(chOk)); pk = scl; kur = scl;
    for iCh = 1 : numel(chOk)
        fb = filterLFP(sig(:, iCh), 'fs', fs, 'type', 'butter', ...
            'dataOnly', true, 'order', 5, 'passband', PASSBAND, ...
            'graphics', false);
        scl(iCh) = 1.4826 * median(abs(fb - median(fb)));
        pk(iCh)  = prctile(abs(fb), 99.9);
        kur(iCh) = kurtosis(fb);
    end
    score = pk ./ scl;

    if isfield(res(iB), 'ED') && ~isempty(res(iB).ED)
        prof = res(iB).ED.prof;
    else
        prof = res(iB).gate.prof;
    end
    edAmp = max(abs(prof), [], 1);

    % rank each criterion's pick within the edAmp ordering (1 = best)
    [~, ordED] = sort(edAmp, 'descend');
    rnk = @(c) find(ordED == find(c == max(c), 1), 1);

    fprintf('\n===== %s =====  (%d ch)\n', basename, numel(chOk));
    fprintf('  p99.9      -> ch %2d, edAmp rank %2d (%.0f of %.0f uV)\n', ...
        chOk(pk == max(pk)), rnk(pk), edAmp(pk == max(pk)), max(edAmp));
    fprintf('  p99.9/scl  -> ch %2d, edAmp rank %2d (%.0f of %.0f uV)\n', ...
        chOk(score == max(score)), rnk(score), ...
        edAmp(score == max(score)), max(edAmp));
    fprintf('  kurtosis   -> ch %2d, edAmp rank %2d (%.0f of %.0f uV)\n', ...
        chOk(kur == max(kur)), rnk(kur), edAmp(kur == max(kur)), max(edAmp));
end
