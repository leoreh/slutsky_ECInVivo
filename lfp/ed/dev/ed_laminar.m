% ED_LAMINAR  Is spatial extent a shape-free criterion for a discharge?
%
% Maslarova et al. 2025 (Nat Commun) report that in mice an IED is visible
% SIMULTANEOUSLY ON EVERY hippocampal channel and layer, whereas a sharp
% wave-ripple is confined to CA1 and its morphology on any one channel is
% laminar - a sharp negative spike in the dendritic layers, a large positive
% slow wave in the pyramidal layer. If that holds here, "how many channels
% carry this event" is a better ED criterion than any waveform template,
% because it does not care which layer the detection channel happens to sit in.
%
% This script measures, for each raMCU mouse and each event class:
%   nChHit - how many probe channels show a deflection > THR robust scales
%   spread - median-over-max of the per-channel peak. AMPLITUDE-INVARIANT: a
%            local event is large on one channel and small elsewhere (-> 0), a
%            globally propagated one is comparable everywhere (-> 1). This is
%            the measure that matters; nChHit alone only says "big".
%   prof   - the per-channel mean waveform (the laminar profile)
%
% Event classes:
%   'ED'    curated discharges          (raMCU3/4/5, from ed_curatedTimes.mat)
%   'gate'  events the current gate accepts (all five mice)
%   'rej'   the STRONGEST gate-rejected candidates - the "just LFP deflections"
%           the curator threw away. The control that matters.
%   'ripp'  accepted ripples            (layer-confined control)
%   'rand'  random times                (null)
%
% Writes ed_laminar.mat next to this file. Run from anywhere.

clear
DEVDIR = fileparts(mfilename('fullpath'));
HALFWIN = 0.15;             % half-window read around each event [s]
THR     = 6;                % channel counts as hit above this many scales
NRAND   = 400;              % random snippets for the per-channel scale
NMAX    = 120;              % cap events per class (read cost)

bps = { ...
    'D:\Data\RA\raMCU1\raMCU1_080621_0930', ...
    'D:\Data\RA\raMCU2\raMCU2_080621_0930', ...
    'D:\Data\RA\raMCU3\raMCU3_211203_084720', ...
    'D:\Data\RA\raMCU4\raMCU4_211220_0834', ...
    'D:\Data\RA\raMCU5\raMCU5_220322_1906'};

cur = load(fullfile(DEVDIR, 'ed_curatedTimes.mat'));
cur = cur.edCurated;
met = ed_methods();

res = struct([]);
for iB = 1 : numel(bps)

    basepath = bps{iB};
    [~, basename] = fileparts(basepath);
    fprintf('\n===== %s =====\n', basename);

    v = basepaths2vars('basepaths', {basepath}, 'vars', {'session'});
    ex = v.session.extracellular;
    nCh = ex.nChannels;
    fs  = ex.srLfp;
    chNeural = sort(unique([ex.spikeGroups.channels{:}]));
    fname = fullfile(basepath, [basename '.lfp']);
    d = dir(fname);
    durRec = d.bytes / 2 / nCh / fs;
    fprintf('  %d ch (%d neural), %.1f h\n', nCh, numel(chNeural), ...
        durRec / 3600);

    % --------------------------------------------------------------- events
    evT = struct();

    iCur = find(strcmp({cur.basename}, basename));
    if ~isempty(iCur)
        evT.ED = cur(iCur).peakTime(:);
    end

    ed = ed_detect(basepath, 'met', met);
    ed.accepted = evt_gate(ed, met.qa);
    evT.gate = ed.peakTime(ed.accepted);
    fprintf('  detect %d cand -> %d gated\n', numel(ed.peakTime), ...
        numel(evT.gate));

    % strongest rejects: the deflections the gate (and the curator) discard
    iRej = find(~ed.accepted);
    [~, iSrt] = sort(ed.fastZ(iRej), 'descend');
    evT.rej = ed.peakTime(iRej(iSrt(1 : min(NMAX, numel(iSrt)))));

    fRipp = fullfile(basepath, [basename '.ripp.mat']);
    if isfile(fRipp)
        r = load(fRipp, 'ripp');
        rt = r.ripp.peakTime(:);
        if isfield(r.ripp, 'accepted') && ~isempty(r.ripp.accepted)
            rt = rt(logical(r.ripp.accepted));
        end
        evT.ripp = rt;
    end

    evT.rand = HALFWIN + rand(NMAX, 1) * (durRec - 2 * HALFWIN);

    % ------------------------------------------------ per-channel robust scale
    tRand = HALFWIN + rand(NRAND, 1) * (durRec - 2 * HALFWIN);
    smp = nan(NRAND, numel(chNeural));
    for iR = 1 : NRAND
        seg = readSeg(fname, tRand(iR), HALFWIN, fs, nCh, chNeural);
        smp(iR, :) = seg(round(size(seg, 1) / 2), :);
    end
    sclCh = 1.4826 * median(abs(smp - median(smp, 1, 'omitnan')), 1, ...
        'omitnan');

    % ------------------------------------------------------------- measure
    cls = fieldnames(evT);
    for iC = 1 : numel(cls)
        t = evT.(cls{iC});
        t = t(t > HALFWIN & t < durRec - HALFWIN);
        if numel(t) > NMAX
            t = t(round(linspace(1, numel(t), NMAX)));
        end
        if isempty(t), continue, end

        nEv = numel(t);
        nChHit = nan(nEv, 1);
        spread = nan(nEv, 1);
        ampMax = nan(nEv, 1);
        prof = nan(round(2 * HALFWIN * fs), numel(chNeural), nEv);
        for iE = 1 : nEv
            seg = readSeg(fname, t(iE), HALFWIN, fs, nCh, chNeural);
            nS = min(size(seg, 1), size(prof, 1));
            base = median(seg(1 : round(0.05 * fs), :), 1);
            seg = seg - base;
            prof(1 : nS, :, iE) = seg(1 : nS, :);
            pk = max(abs(seg), [], 1);
            nChHit(iE) = sum(pk(:)' > THR * sclCh);
            pkZ = pk(:)' ./ sclCh;          % per-channel, in its own units
            spread(iE) = median(pkZ) / max(pkZ);
            ampMax(iE) = max(pkZ);
        end

        res(iB).(cls{iC}).nChHit = nChHit;
        res(iB).(cls{iC}).spread = spread;
        res(iB).(cls{iC}).ampMax = ampMax;
        res(iB).(cls{iC}).prof   = mean(prof, 3, 'omitnan');
        res(iB).(cls{iC}).t      = t;
        fprintf(['  %-5s n=%3d  nChHit %4.1f/%d   spread %.2f   ', ...
            'ampMax %5.1f\n'], cls{iC}, nEv, median(nChHit), ...
            numel(chNeural), median(spread), median(ampMax));
    end

    res(iB).basename = basename;
    res(iB).chNeural = chNeural;
    res(iB).sclCh    = sclCh;
    res(iB).fs       = fs;
end

save(fullfile(DEVDIR, 'ed_laminar.mat'), 'res', '-v7.3');
fprintf('\nsaved ed_laminar.mat\n');


% =========================================================================
%  LOCAL
% =========================================================================
function seg = readSeg(fname, t, halfWin, fs, nCh, ch)
seg = binary_load(fname, 'start', t - halfWin, 'duration', 2 * halfWin, ...
    'fs', fs, 'nCh', nCh, 'ch', ch, 'bit2uv', 0.195);
seg = double(seg);
end
