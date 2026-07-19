% EVT_DOC  Event pipelines (ripples + ED) — shared architecture and file map.
%
% Two detection pipelines — sharp-wave ripples (lfp/ripples) and epileptiform
% discharges (lfp/ed) — run in parallel on a shared, event-agnostic layer
% (lfp/events, the evt_* functions) and hand a common curation-ready struct to
% the GUI (guiPath). One session runs end to end through a wrapper; every
% result is a <basename>.<var>.mat file. This file is the shared-design
% reference for both pipelines.
%
% The two pipelines are parallel by construction: they share every
% detection-agnostic step (spike prep, bout times, QA, matched controls, state
% labelling, maps, spike analysis, saving, plotting) and differ only where the
% modality genuinely differs — the detection model, the per-event features, and
% the signal loaded. Only detection is modality-specific: ripples z-score a
% narrowband envelope; ED z-score a moving baseline to flag sharp transients.
% Every side effect is gated (flgSave writes disk, flgPlot draws figures,
% flgCurate opens the GUI), so a pipeline runs headless for verification.
%
% # Entry points
%
% ripp_wrapper  Ripple orchestrator. Reads top-to-bottom: setup -> detect-or-
%               load -> curate. Detects on one hippocampal LFP channel.
% ed_wrapper    ED orchestrator, same skeleton. Detects EEG or LFP.
%
% Re-run model (both): if <basename>.<name>.mat exists and flgForce is false,
% the stored result is loaded and handed to the GUI (the cheap re-curate path);
% otherwise the full pipeline runs. Signal loading and the analyses live inside
% the detect branch, so re-curation touches no signals.
%
% # Canonical event schema
%
% Both ripp (<basename>.ripp.mat) and ed (<basename>.ed.mat) carry the same
% curation contract plus per-modality metric columns:
%   .times     [N x 2]  event start / stop (s, absolute)
%   .peakTime  [N x 1]  event peak (s, absolute)
%   .state     [N x 1]  categorical vigilance state (<undefined> if none)
%   .accepted  [N x 1]  logical curation mask (seeded all-true post-QA)
%   .ctrlTimes [N x 2]  matched control intervals (analysis)
%   .info      struct   detection params, fs, basename, win, runtime
% Ripple metric columns: .amp .freq .freqEvent .energy .dur .skew .emg .spkGain.
% ED metric columns:     .amp .ampZ .dur .width10 .emgZ.
% Automatic QA runs as a filter at detection (evt_qa + evt_subset): only
% survivors are saved and .accepted starts all-true; curation flips it.
%
% # Shared layer (lfp/events)
%
% Event-agnostic; both wrappers call these with modality-specific inputs.
% evt_files     Standard output file paths for a modality ('ripp' / 'ed').
% evt_spkPrep   Window-relative single-unit + pooled-MUA times, unit types.
% evt_boutTimes Window-relative bout times + valid-state + NREM baseline sets.
% evt_emgScore  Per-event EMG z vs a baseline window (a QA metric).
% evt_spkGain   Per-event MUA spike-gain z (a QA metric; ripples).
% evt_qa        Combine QA criteria (state inclusion + metric ranges) into one
%               pass mask; the wrapper drops the failures via evt_subset.
% evt_subset    Subset every per-event field by a logical mask.
% evt_ctrlTimes Duration-matched control intervals from valid states.
% evt_states    Per-event vigilance state + the per-bout rate/density table.
% evt_maps      Per-event signal maps around each event peak.
% evt_spks      Spike entry point: per-unit stats + per-event population metrics
%               + 3D raster + per-unit PETH (orchestrates evt_spksParams,
%               evt_spkPeth, evt_pethNorm, evt_rankOrder).
% evt_saveSpks  Split the spike result into a light stats file + heavy raster.
% evt_plotSpks  Spike-modulation summary figure.
% evt2ns        Export events to a NeuroScope .evt file (start / peak / stop
%               marks, labelled by acceptance). Ripples wire it via flgNS.
% guiPath Manual curation; reads and writes the .accepted mask.
%
% # Modality-specific (deliberately not shared)
%
% Ripples: ripp_pickCh (detection channel), ripp_sigLoad (channel + EMG),
%   ripp_sigPrep (filter + envelope + z-score), ripp_times (threshold
%   candidates), ripp_params (frequency / energy / skew), ripp_detect (the
%   shared detect -> params -> maps -> QA core), spklfp_phase (spike-LFP).
% QA marks .accepted (ripples, never removed) or subsets (ED); mcu_tblVivo
%   filters .accepted for the per-event ripple table.
% ED: ed_sigLoad (EEG or LFP channel + EMG), ed_detect (moving-z transients),
%   ed_params (half- / 10%-amplitude widths).
% The channel loaders share ripp_pickCh (the ripple channel) and evt_loadCh
% (binary .lfp read + averaging + bit2uv autodetect).
%
% # Outputs  (<basename>.<var>.mat, written when flgSave = true)
%
%   Concern             Ripples             ED                Producer
%   Events (curation)   .ripp.mat           .ed.mat           wrapper
%   Per-bout states     .rippStates.mat     .edStates.mat     evt_states
%   Signal maps         .rippMaps.mat       .edMaps.mat       evt_maps
%   Spikes: stats+PETH  .rippSpks.mat       .edSpks.mat       evt_saveSpks
%   Spike raster (3D)   .rippSpkMaps.mat    .edSpkMaps.mat    evt_saveSpks
%   Phase coupling      .rippSpkLfp.mat     —                 spklfp_phase
%
% The light spikes file (.rippSpks / .edSpks) carries per-unit scalar stats +
% the per-unit PETH + tstamps; the manuscript loads it hot. The 3D raster is
% split into its own heavy file. A per-type population PETH is not precomputed;
% derive it from the 3D raster if a table ever needs it.
%
% # Robustness
%
% Both wrappers degrade gracefully when sleep states, spikes, units, or EMG are
% absent (real EA / RA sessions): missing spikes skip the spike analyses;
% missing EMG / NREM relax the corresponding QA criteria; 'nrem' z-scoring falls
% back to 'adaptive'.
%
% Terminology: a "map" or "raster" is the 3D per-event array [unit x event x
% bin]; a "PETH" is its per-unit average over events (2D [unit x bin]).
