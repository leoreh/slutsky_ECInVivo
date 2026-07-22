% MCU_ED  Epileptiform discharges: detect -> curate by cluster -> count.
%
% A scratch pad, not a function. Run a section at a time. The pipeline itself is
% in lfp/ed; see lfp/ed/dev/ed_pipeline_rebuild.md for how every default was
% measured, and evt_doc for how the ED and ripple pipelines share their spine.
%
% The shape mirrors mcu_ripples: loop 1 detects every session unattended, loop 2
% is manual and one mouse at a time, and the last sections read the result.
%
% What curation is now. Detection is permissive and proposes thousands of
% candidates on a 24 h recording. met.qa cuts that to a POOL of a few hundred by
% asking only what no discharge can fail - is it sharp, does it stand alone -
% and ed_curate groups the pool by waveform shape so you accept whole TYPES.
% Nothing in the pipeline encodes what a discharge looks like; you do, per
% mouse, by ticking clusters.
%
% Requires a binary <basename>.lfp and a session.mat: detection reads one
% auto-picked raw channel (ed_pickCh). The EA cohort has neither and is no
% longer runnable here.
%
% CAUTION: raMCU3 / raMCU4 / raMCU5 hold hand-curated masks. Re-detecting them
% overwrites the mask (the old file goes to bkup/ first, and the curated peak
% times are also kept in lfp/ed/dev/ed_curatedTimes.mat).



%% ========================================================================
%  DETECT  (loop 1: every session, unattended)
%  ========================================================================
% ~20 s per 24 h session, most of it the channel probe. Writes
% <basename>.ed.mat (all candidates + metrics), .edMaps.mat (one waveform per
% candidate) and, via the headless filter, .edStates.mat. flgForce backs up any
% existing file first, so a mask you already curated is recoverable from bkup/.

basepaths = [mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl'), ...
    mcu_basepaths('ra')];
nFiles = numel(basepaths);
met = ed_methods('default');

for iFile = 1 : nFiles
    ed_wrapper('basepath', basepaths{iFile}, 'met', met, 'win', [0 Inf], ...
        'flgSave', true, 'flgForce', true);
end


%% ========================================================================
%  CURATE  (loop 2: manual, one mouse at a time)
%  ========================================================================
% Tick the clusters whose median waveform is a discharge, then Save. Reading
% the tiles:
%   - a discharge is a sharp complex ~20-30 ms wide on flat background, back to
%     baseline within 50-100 ms. It may be positive-going or negative-going;
%     polarity depends on which layer the electrode sits in, so do NOT reject a
%     cluster for having the "wrong" sign.
%   - an ordinary sharp wave is a slow monophasic negative excursion peaking
%     ~40 ms out and decaying over 300+ ms. This is the bulk of the pool.
%   - a step artifact drops and never returns.
%   - the MUA row is a second opinion, not a criterion: a discharge cluster
%     usually shows firing falling well below baseline for 50-300 ms after the
%     event, while artifacts and ripples sit flat at 1.
% 'clusters' re-runs with a different count. More clusters means finer types
% and more boxes; 12 measured best, but a mouse with a big pool may want more.

iFile = 6;

ed_curate(basepaths{iFile});

% optional per-event pass over what you accepted, for a mouse you are unsure
% about. Up / Down accept-reject, Left / Right step, Ctrl+S save.
[~, vm, gm] = guiPath(basepaths{iFile}, 'preset', 'ed');

% reopen fast (the signals stay loaded; only the events are re-read)
vm.ed.data = [];
guiPath(basepaths{iFile}, 'varMap', vm, 'guiMap', gm);


%% ========================================================================
%  RATE BY STATE  (the reported number)
%  ========================================================================
% One row per session x state, so ONE NUMBER PER MOUSE PER STATE: accepted
% discharges divided by that state's scored exposure, in events per minute.
% The mouse is the unit genotype was assigned at, so it is the unit that
% enters the figure and the model - pooling events across mice would let a
% long recording or a busy animal count as several.
%
% 'ALL' is every accepted event over the whole recording and is NOT the sum of
% the state rows: an event whose peak falls in an unscored gap belongs to no
% state, so the gap between them reads how completely the session was scored.
%
% States are merged for reporting (evt_stateMerge): QWAKE into WAKE, LSLEEP
% into NREM. The scoring on disk is untouched. In the lh cohort those two carry
% 38-121 min each, so the merge is not cosmetic - it moves WAKE and NREM
% exposure by roughly 20%. N/REM survives as its own row and usually falls
% below the exposure floor below; merge it too if you would rather it did not
% appear at all.

tblEd = ed_tbl(basepaths);
tblEd.genotype = mcu_geno(tblEd.sbjID);

% whole-session burden per mouse
tblAll = tblEd(tblEd.state == 'ALL', :);
guiTbl_bar(tblAll, 'xVar', 'genotype', 'yVar', 'edRate');

% state dependence. Drop ALL, and drop states with too little exposure to give
% a stable rate - a REM bout total of two minutes turns one event into 0.5/min.
% lme_analyse drops empty categorical levels itself, so no removecats here.
tblState = tblEd(tblEd.state ~= 'ALL' & tblEd.durState > 30, :);
guiTbl_bar(tblState, 'xVar', 'state', 'yVar', 'edRate', 'grpVar', 'genotype');

% -> PRISM, grouped layout: row per state, column per genotype, one subcolumn
% per mouse. Padded to the widest genotype, so blanks are missing mice.
%
% flgSort false keeps the CATEGORY order instead of sorting alphabetically,
% which is the only way the columns come out Control / MCU-KO / CAG-MCU-KO
% rather than CAG first. Same for the state rows, hence the reorder.
tblState.state = removecats(tblState.state);
ssCfg = as_loadConfig([]);              % WAKE QWAKE LSLEEP NREM REM ...
tblState.state = reordercats(tblState.state, ...
    intersect(ssCfg.names, categories(tblState.state), 'stable'));

tbl2prism(tblState, 'yVar', 'edRate', 'grpVar', 'genotype', ...
    'rowVar', 'state', 'flgSort', false);

% counts per mouse, same layout - the denominator-free number to report
% alongside the rate. Copy one, paste, then run the other.
tbl2prism(tblState, 'yVar', 'nEd', 'grpVar', 'genotype', ...
    'rowVar', 'state', 'flgSort', false);

% CHECK BEFORE MODELLING. If the discharges are confined to the CAG cohort,
% every Control and MCU-KO rate is 0 and the interaction below is fit against
% two columns of zeros: the p-values are not wrong so much as meaningless, and
% a presence/absence test (any discharge at all, Fisher) is the honest claim.
% This prints, per genotype, how many mice carry any discharge at all.
disp(varfun(@(x) [nnz(x > 0), numel(x)], tblAll, ...
    'GroupingVariables', 'genotype', 'InputVariables', 'edRate'));

% Modelled as a rate, with a floor on exposure. The statistically cleaner form
% is a count model with log(durState) as an OFFSET - a state carrying 20 min of
% REM and one carrying 8 h of NREM should not weigh the same - but lme_fit has
% no offset argument, so that would need wiring in. The durState filter above is
% the blunt stand-in; check it before reading a REM effect.
frml = 'edRate ~ state * genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblState, frml);


%% ========================================================================
%  WAVEFORM PER MOUSE
%  ========================================================================
% One row per ACCEPTED event, carrying its detrended snippet. Tile by genotype
% and group by mouse (or the reverse) and the legend gives the event count per
% group - which is the per-mouse n, read off the same figure that shows whether
% the waveforms agree.
%
% Detrended, not normalised: amplitude is real here and worth seeing. Set
% Dispersion to Spread and Stat to Median for a robust central trace.
%
% Aligned per session (ed_wvTbl 'align', default 'auto'). Detection centres
% each event on its largest absolute swing, so some sit on the peak and some on
% the trough and the average smears; aligning every event to one feature
% sharpens t=0. 'auto' picks that feature from each mouse's OWN average - trough
% for a negative mouse, peak for a positive one - since polarity depends on the
% recording layer. Force it with 'align','trough'|'peak', or 'none' to see the
% detection alignment.

[tblWv, tstamps] = ed_wvTbl(basepaths, 'align', 'trough');
tblWv.genotype = mcu_geno(tblWv.sbjID);

guiTbl_xy(tstamps * 1000, tblWv, 'yVar', 'lfp', 'tileVar', 'genotype', ...
    'grpVar', 'sbjID', 'xLbl', 'time (ms)', 'xLim', [-50 50]);

% counts per mouse as a table, if the legend is not enough
tblN = groupsummary(tblWv, {'genotype', 'sbjID'});

% -> PRISM, ONE XY block per genotype: X = time, one column per mouse carrying
% Mean / SD / N (across that mouse's events) at each time point. In Prism make
% an XY table, Format Data Table -> "Enter and plot error... Mean, SD, N", and
% paste - the pasted block's first row is the mouse names, one over each triple.
% wv2prism copies the level named in 'copy'; re-run with the next to paste it
% into its own graph. Same detrended, trough-aligned waveforms as the GUI above.
wvBlocks = wv2prism(tblWv, tstamps * 1000, 'grpVar', 'sbjID', ...
    'splitVar', 'genotype', 'xLbl', 'time (ms)', 'xLim', [-50 50], ...
    'copy', 'CAG-MCU-KO');



