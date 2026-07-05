% ED_DOC  Epileptiform-discharge (ED) pipeline — architecture and file map.
%
% The pipeline detects epileptiform discharges (EDs) — sharp interictal spikes
% — on one EEG or LFP channel, characterizes each event, relates events to
% single- and multi-unit spiking and to vigilance state, and hands a
% curation-ready struct to the shared events GUI. One session runs end to end
% through ed_wrapper; every result is a <basename>.<var>.mat file. The pipeline
% mirrors the ripple pipeline (see lfp/ripples/ripp_doc) and shares its
% event-agnostic layer.
%
% Implementation is MATLAB on the shared event layer (lfp/events, the evt_*
% functions), so the ED and ripple pipelines stay parallel and share the
% detection-agnostic steps: matched controls, state labelling, maps, spike
% analysis, and plotting. Only detection is ED-specific — a moving-baseline
% z-score flags sharp transients (positive, negative, or both) rather than a
% narrowband envelope. Every side effect is gated (flgSave writes disk, flgPlot
% draws figures, flgCurate opens the GUI), so the pipeline runs headless for
% verification. If <basename>.ed.mat exists and flgForce is false, detection is
% skipped and the stored result loads straight into the GUI (the re-curate path).
% The shared layer is documented in lfp/events/REFRACTOR.md; this file covers
% the ED-specific scripts.
%
% # Scripts
%
% ed_wrapper    Orchestrator and only entry point: setup -> signal -> detect ->
%               characterize -> QA filter -> parity analyses -> save / plot /
%               curate.
% ed_sigLoad    Loads the detection signal (EEG sSig.eeg, or a raw LFP channel),
%               the EMG (z-score and RMS forms), a spectrogram adapter, and the
%               full-session sSig for the curation GUI.
% ed_detect     Detection. Flags sharp transients where the signal crosses a
%               moving-baseline z-score threshold, merging twin peaks and
%               enforcing a refractory window.
% ed_params     Per-event features: amplitude (.amp), z-scored amplitude
%               (.ampZ), half-amplitude width (.dur), and 10%-width (.width10).
% ed_reject_emg QA scoring: flags high-EMG events by z-score or emg_rms, writing
%               the .idxQA.emg mask that ed_wrapper filters on.
%
% Shared event layer (lfp/events), called by ed_wrapper:
% evt_spkPrep   Window-relative single-unit + pooled-MUA spike times, unit types.
% evt_ctrlTimes Duration-matched control intervals drawn from valid states.
% evt_states    Per-event vigilance state, plus the per-bout rate/density table.
% evt_maps      Per-event signal maps around each discharge peak.
% evt_spks      Spike entry point: per-unit stats + per-event population metrics
%               + 3D raster + per-unit PETH (orchestrates evt_spksParams,
%               evt_spkPeth, evt_pethNorm).
% evt_plotSpks  Spike-modulation summary figure.
% evt_viewSpks  Interactive map / PETH viewers; population PETH on demand
%               (evt_pethPop) from the 3D raster.
% evt_rate      Optional discharge-rate time series (flgRate).
% gui_curate    Manual curation; reads and writes the .ed.mat .accepted mask.
%
% # Outputs  (<basename>.<var>.mat, written when flgSave = true)
%
% ed          Events + per-event metrics: .times .peakTime .pos .state
%             .accepted .ctrlTimes, discharge params (.amp .ampZ .dur .width10),
%             QA metric (.emg), and population spike metrics (.spks).
% edStates    Per-bout table — discharge Rate / Density / Duration by vigilance
%             state (when sleep states are present).
% edMaps      Per-event signal maps (.lfp) + .tstamps.
% edSpks      Per-unit spike stats (frEvt, frCtrl, frZ, com, asym, rankMean, …)
%             + the per-unit PETH (.peth) + .tstamps. Light.
% edSpkMaps   The 3D spike raster [unit x event x bin] (.su / .mu, each
%             .evt / .ctrl) + .tstamps. Heavy — the source for any PETH reduction.
%
% Terminology: a "map" or "raster" is the 3D per-event array [unit x event x
% bin]; a "PETH" (peri-event time histogram) is its per-unit average over
% events (2D [unit x bin]).
