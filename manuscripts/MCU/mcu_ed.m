% MCU_ED  Epileptiform discharges: detect -> curate -> count by state.
%
% A scratch pad, not a function. Run a section at a time. The pipeline itself is
% in lfp/ed; see lfp/ed/dev/ed_pipeline_rebuild.md for how every default was
% measured, and evt_doc for how the ED and ripple pipelines share their spine.
%
% The shape mirrors mcu_ripples: loop 1 detects every session, loop 2 is manual
% and one mouse at a time, and the last sections read the result. Detection is
% deliberately permissive (thousands of candidates on a 24 h recording), so the
% curation pass is not optional - it is where the operating point is chosen.



%% ========================================================================
%  DETECT  (loop 1: every session, unattended)
%  ========================================================================
% ~12 s per 24 h session. Writes <basename>.ed.mat (all candidates + metrics),
% .edMaps.mat (one waveform per candidate) and, via the headless gate,
% .edStates.mat. flgForce backs up any existing file before overwriting, so a
% mask you already curated is recoverable from bkup/.

basepaths = [mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl'), ...
    mcu_basepaths('ra')];
nFiles = numel(basepaths);
met = ed_methods('default');        % detection + default QA filter (met.qa)

for iFile = 12 : nFiles
    ed_wrapper('basepath', basepaths{iFile}, 'met', met, 'win', [0 Inf], ...
        'flgSave', true, 'flgForce', true);
end


%% ========================================================================
%  CURATE  (loop 2: manual, one mouse at a time)
%  ========================================================================
% Two passes. The first sets four thresholds over ALL candidates at once
% and shows the kept-vs-removed mean waveform per state; the second steps
% the survivors one by one. Doing the second without the first means
% walking thousands of events.

iFile = 15;

% pass 1 - bulk. Move a threshold, watch the waveform split, press Save.
ed_curate(basepaths{iFile}, 'qa', met.qa);

% pass 2 - per event. Up / Down accept-reject, Left / Right step, Ctrl+S save.
[~, vm, gm] = guiPath(basepaths{iFile}, 'preset', 'ed');

% reopen fast (the signals stay loaded; only the events are re-read)
vm.ed.data = [];
guiPath(basepaths{iFile}, 'varMap', vm, 'guiMap', gm);


%% ========================================================================
%  RE-GATE  (optional: apply one spec to every session, no re-detection)
%  ========================================================================
% Curation is separable from detection, so a threshold can be revisited across
% the cohort in seconds. This OVERWRITES the per-mouse masks set in the GUI -
% run it before the manual pass, not after.

qa = met.qa;
qa.ranges.ampG = [10 Inf];          % e.g. tighten the amplitude criterion

for iFile = 1 : nFiles
    ed_curate(basepaths{iFile}, 'qa', qa, 'flgGui', false);
end


%% ========================================================================
%  COUNTS BY STATE
%  ========================================================================
% One row per session x state. 'ALL' is every accepted event over the whole
% recording and is NOT the sum of the state rows - an event whose peak falls in
% an unscored gap belongs to no state, so the gap between them reads how
% completely the session was scored.

tblEd = ed_tbl(basepaths);

tblEd.genotype = mcu_geno(tblEd.sbjID);

% whole-session burden per mouse
tblAll = tblEd(tblEd.state == 'ALL', :);
guiTbl_bar(tblAll, 'xVar', 'genotype', 'yVar', 'edRate');

% state dependence (drop ALL, and states with little exposure). lme_analyse
% drops empty categorical levels itself, so no removecats is needed here.
tblState = tblEd(tblEd.state ~= 'ALL' & tblEd.durState > 0.5, :);
guiTbl_bar(tblState, 'xVar', 'state', 'yVar', 'edRate', 'grpVar', 'genotype');

% Modelled as a rate, with a floor on exposure. The statistically cleaner form
% is a count model with log(durState) as an OFFSET - a state carrying 20 min of
% REM and one carrying 8 h of NREM should not weigh the same - but lme_fit has
% no offset argument, so that would need wiring in. The durState filter above is
% the blunt stand-in; check it before reading a REM effect.
frml = 'edRate ~ state * genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblState, frml);


%% ========================================================================
%  SANITY: WHAT DID THE GATE KEEP?
%  ========================================================================
% Worth one look per cohort. A non-epileptic mouse should land near a few tens
% of events; lh132 currently returns an epileptic-range rate, which is either
% real or a bad channel.

iFile = 1;
[~, basename] = fileparts(basepaths{iFile});
load(fullfile(basepaths{iFile}, [basename, '.ed.mat']), 'ed');

acc = ed.accepted;
fprintf('%s: %d / %d accepted\n', basename, nnz(acc), numel(acc));
fprintf('  ampG %.1f | ampZ %.1f | hfRatio %.2f | dur %.1f ms | emg %.2f\n', ...
    median(ed.ampG(acc), 'omitnan'), median(ed.ampZ(acc), 'omitnan'), ...
    median(ed.hfRatio(acc), 'omitnan'), median(ed.dur(acc), 'omitnan'), ...
    median(ed.emg(acc), 'omitnan'));

% the four metrics against each other, kept vs removed
tblEvt = table(ed.ampG, ed.ampZ, ed.hfRatio, ed.emg, ed.dur, ...
    categorical(acc, [false true], {'removed', 'kept'}), ...
    'VariableNames', {'ampG', 'ampZ', 'hfRatio', 'emg', 'dur', 'status'});
guiTbl_scatHist(tblEvt, 'xVar', 'ampG', 'yVar', 'hfRatio', 'grpVar', 'status');


%% ========================================================================
%  THE EA COHORT  (epileptic mice; flat folder, no session.mat / sleep_states)
%  ========================================================================
% Several recordings share one directory and the stem does not follow it, so
% every entry point takes an explicit basename. With no sleep scoring the events
% are left unlabelled, ed_tbl returns only the 'ALL' row, and guiPath drops the
% state strip - the rest of the pipeline is unchanged.

eaPath  = 'D:\Data\EA';
eaNames = {'220611_0750', '220615_0801', '220824_0906'};

for iFile = 1 : numel(eaNames)
    ed_wrapper('basepath', eaPath, 'basename', eaNames{iFile}, 'met', met, ...
        'flgSave', true, 'flgForce', true);
end

ed_curate(eaPath, 'basename', eaNames{1}, 'qa', met.qa);
guiPath(eaPath, 'basename', eaNames{1}, 'preset', 'ed');

tblEA = ed_tbl(repmat({eaPath}, 1, numel(eaNames)), eaNames);

% EOF
