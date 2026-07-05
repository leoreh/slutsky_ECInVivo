% RIPP_DOC  Sharp-wave ripple (SWR) pipeline — architecture and file map.
%
% The pipeline detects sharp-wave ripples (SWRs) on one hippocampal LFP
% channel, characterizes each event, relates events to single- and multi-unit
% spiking and to vigilance state, and hands a curation-ready struct to the
% shared events GUI. One session runs end to end through ripp_wrapper; every
% result is a <basename>.<var>.mat file.
%
% Implementation is MATLAB built on the shared event layer (lfp/events, the
% evt_* functions), so the ripple and epileptiform-discharge (ED) pipelines
% stay parallel and share the detection-agnostic steps: matched controls,
% state labelling, maps, spike analysis, and plotting. Only detection is
% ripple-specific — a narrowband envelope is z-scored and thresholded. Every
% side effect is gated (flgSave writes disk, flgPlot draws figures, flgCurate
% opens the GUI), so the pipeline runs headless for verification. The shared
% layer is documented in lfp/events/REFRACTOR.md; this file covers the
% ripple-specific scripts.
%
% # Scripts
%
% ripp_wrapper  Orchestrator and only entry point. Reads top-to-bottom as the
%               pipeline: setup -> signal -> detect -> characterize -> spikes
%               -> phase -> save / plot / curate.
% ripp_sigPrep  Bandpass-filters the LFP, takes the analytic (Hilbert)
%               envelope, and z-scores it against an NREM or adaptive baseline.
%               Also called by the curation GUI to redraw ripple signals.
% ripp_times    Thresholds the z-scored signal (start / peak / continuation
%               thresholds, duration limits, event merging) into candidates.
% ripp_qa       Returns a pass mask from three criteria — vigilance state, MUA
%               spike-gain, and EMG — each skipped when its data is absent.
%               ripp_wrapper drops the failures.
% ripp_params   Per-event features from the analytic signal: duration,
%               amplitude, frequency, energy, skewness.
% ripp2ns       Exports events to Neuroscope .evt files for external inspection.
%
% Shared event layer (lfp/events), called by ripp_wrapper:
% evt_spkPrep   Window-relative single-unit + pooled-MUA spike times, unit types.
% evt_ctrlTimes Duration-matched control intervals drawn from valid states.
% evt_states    Per-event vigilance state, plus the per-bout rate/density table.
% evt_maps      Per-event LFP maps around each ripple peak.
% evt_spks      Spike entry point: per-unit stats + per-event population metrics
%               + 3D raster + per-unit PETH (orchestrates evt_spksParams,
%               evt_spkPeth, evt_pethNorm).
% evt_plotSpks  Spike-modulation summary figure.
% evt_viewSpks  Interactive map / PETH viewers; population PETH computed on
%               demand (evt_pethPop) from the 3D raster.
% spklfp_phase  Spike-LFP phase coupling.
% gui_curate    Manual curation; reads and writes the .ripp.mat .accepted mask.
%
% # Outputs  (<basename>.<var>.mat, written when flgSave = true)
%
% ripp        Events + per-event metrics: .times .peakTime .state .accepted
%             .ctrlTimes, ripple params (.amp .freq .freqEvent .energy .dur
%             .skew), QA metrics (.emg .spkGain), population spike metrics (.spks).
% rippStates  Per-bout table — ripple Rate / Density / Duration by vigilance
%             state (the "SWR rate by brain state" analysis).
% rippMaps    Per-event LFP maps (.lfp .filt .amp .freq .z) + .tstamps.
% rippSpks    Per-unit spike stats (frEvt, frCtrl, frZ, com, asym, rankMean, …)
%             + the per-unit PETH (.peth) + .tstamps. Light — the manuscript
%             loader (mcu_tblVivo) reads this file.
% rippSpkMaps The 3D spike raster [unit x event x bin] (.su / .mu, each
%             .evt / .ctrl) + .tstamps. Heavy — the source for any PETH
%             reduction (per-unit or population).
% rippSpkLfp  Spike-LFP phase coupling.
%
% Terminology: a "map" or "raster" is the 3D per-event array [unit x event x
% bin]; a "PETH" (peri-event time histogram) is its per-unit average over
% events (2D [unit x bin]).
