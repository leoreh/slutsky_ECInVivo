function guiPath_doc()

% Walk through the guiPath curation framework end to end.
%
% This file is documentation you can run: every code line below is a call you
% can paste. guiPath_doc itself does nothing. Read it top to bottom, or type
% help guiPath_doc . Start with CONCEPTS for the ideas, WALKTHROUGH for a live
% session, REFERENCE for the field lists.
%
% The framework shows one session's signals. Several event / state sets can be
% loaded at once; the CURATE selector picks which one is the editable target -
% accept / reject each event, or a state label per epoch - or None, to only
% view. PRESET picks the arrangement (which panels), independently of CURATE.
% Everything rests on one split. WHAT to show is a declarative description you
% build and edit; HOW to fetch and draw it is the framework's job. You describe
% panels and addresses; guiPath_load and guiPath do the loading and drawing.
%
%
% CONCEPTS
%
%   Two configs describe a session, held and passed separately.
%   - cfgData   the panels: what is shown, and the one target you curate. A
%               flat struct, one field per panel. This is what you edit.
%   - cfgGui    the behaviour: display name, session file, mode, window width,
%               save target. Everything guiPath needs beyond the panels.
%   guiPath_presets(name) returns both. guiPath_load fills data into cfgData.
%   guiPath draws them.
%
%   The panel is the atomic unit, one self-describing struct from guiPath_panel.
%   It states WHAT it is (type), WHERE it sits (region), and WHERE its data
%   lives (src). A cfgData is a struct of these, one per field; the field name
%   is the panel name and the field order is the stacking order in the region.
%
%   The address (src) is an indirection. A panel does not hold its data or know
%   how to read it; it holds a short string saying where the data lives, such
%   as 'sleep_sig:eeg' or 'ed'. guiPath_src is the single place that knows how to
%   turn an address into raw data. So a panel stays a plain description, and
%   every way of fetching data lives in one function (see the grammar in
%   REFERENCE).
%
%   Slim versus full. A cfgData from guiPath_presets is slim: addresses only, no
%   data. guiPath_load reads each address and writes the result back into the
%   panel (.data, .fs, a resolved .ylim), returning the same struct now full. A
%   panel that already carries data is skipped, so topping up a config is cheap.
%
%   CURATE picks the editable target. Several event / state sets can be loaded
%   (a preset brings one; Load... adds more) and all are drawn at once; the
%   CURATE selector chooses which single one you edit, and its type sets the
%   mode:
%   - an eventTicks set gives events mode: accept or reject each event.
%   - a stateStrip set gives states mode: assign a state to each epoch.
%   - None gives view mode: no target, Prev / Next just steps the window.
%   A hypnogram and a stateStrip both draw a coloured state strip, but a
%   hypnogram is always read-only context (committed bouts); a stateStrip can be
%   the target.
%
%
% CONCEPTS: ctx, the session cache
%
%   ctx is a small struct that travels with a load. guiPath_ctx(basepath,
%   basename) builds it. It carries the session location (.basepath, .basename)
%   and one cache (.cache) shared by every address resolved during the load.
%
%   The problem it solves. Many panels read from the same file. In the EDs
%   preset the hypnogram, the spectrogram, the EMG, and the EMG RMS all come out
%   of the assembled sleep signals, which ed_sigLoad builds once from the raw
%   LFP (an expensive read). Resolving each address on its own would rebuild
%   that bundle several times. The cache stores each file (and each computed
%   source) the first time it is read and returns the stored copy on every later
%   hit, so a given file is read at most once per ctx.
%
%   The mechanism (worth understanding, because MATLAB makes it surprising). A
%   struct is a value: passing it into a function copies it, and writes inside
%   that function never reach the caller. If .cache were an ordinary struct or
%   array, a cache write inside guiPath_src would be lost on return.
%   containers.Map is one of MATLAB's few reference (handle) types. Copying the
%   ctx struct copies the handle, not the map behind it, so every copy points at
%   the same underlying map. A write through any copy (ctx.cache(key) = value)
%   is seen through all of them. That is why the loaders can take ctx by value,
%   write into ctx.cache, never return ctx, and still leave the caller's cache
%   filled: the struct is copied, the cache is shared.
%
%   Two levels of "load once" stack on top of each other.
%   - Panel level: guiPath_load skips a panel that already has data, and
%     guiPath_presets copies loaded data across matching panels on a preset
%     switch. This avoids re-resolving an address whose result is already held.
%   - File level: even when a panel must load, ctx.cache makes the underlying
%     file read happen once. Two panels with different addresses into the same
%     file ('sleep_sig:emg' and 'sleep_sig:emg_rms') share the one cached read.
%
%   Use case, a preset switch. guiPath builds one ctx when it opens and keeps
%   it for the whole session (in hFig.UserData.ctx). Switching EDs -> Ripples
%   calls guiPath_load again with the SAME ctx. The panels shared by both presets
%   (hypnogram, spectrogram, EMG, EMG RMS, unit raster) are already cached, so
%   only the ripple-specific reads run (the ripp file, the ripple-band LFP).
%   Without a shared ctx each switch would rebuild the sleep signals from raw
%   LFP.
%
%   By hand the cache is usually implicit. guiPath_load(cfgData, basepath) makes
%   a fresh ctx for that one call, so a file is read once within it. Pass your
%   own ctx only when you want to share a cache across several calls:
%       ctx = guiPath_ctx(basepath, basename);
%       cfgA = guiPath_load(cfgA, basepath, ctx);
%       cfgB = guiPath_load(cfgB, basepath, ctx);   % reuses what cfgA read
%
%
% WALKTHROUGH
%
%   Set a session. guiPath takes a basepath and derives basename from it;
%   with no argument it uses the current folder (pwd).
%       basepath = 'D:\Data\lh100\lh100_220413_111004';
%       [~, basename] = fileparts(basepath);
%
%   1. Open a session. A preset is auto-detected from the files present.
%       guiPath(basepath);
%       guiPath(basepath, 'preset', 'States');       % force one by name
%
%   2. Get a preset's panels and behaviour.
%       [cfgData, cfgGui] = guiPath_presets('EDs');   % Ripples States template
%       fieldnames(cfgData)'         % panel names, in stacking order
%       guiPath_presets()             % the preset list {name, file}
%
%   3. Read one panel and its address.
%       cfgData.spec                 % a spec panel in the top region
%       cfgData.spec.type            % 'spec'
%       cfgData.spec.region          % 'top'
%       cfgData.spec.src             % 'fn:spec' (the address)
%
%   4. Load the data. Every address is resolved into .data; a panel that already
%      has data is left alone.
%       cfgData = guiPath_load(cfgData, basepath);
%       cfgData.spec.data            % now present
%
%   5. Open your own config. A slim cfgData (no data) is loaded for you; omit
%      cfgGui and it is derived from the panels.
%       guiPath(basepath, 'cfgData', cfgData, 'cfgGui', cfgGui);
%
%   6. Edit or add a panel. guiPath_panel(type, region, address) builds one; its
%      field position sets its stacking order in the region.
%       cfgData.emg.height = 1.0;
%       cfgData.emg2 = guiPath_panel('trace', 'bottom', 'sleep_sig:emg', ...
%           'label', 'EMG 2');
%       guiPath(basepath, 'cfgData', cfgData);
%
%   7. Build from a blank template.
%       [cfgData, cfgGui] = guiPath_presets('template');  % empty cfgData
%       cfgData.spec = guiPath_panel('spec',  'top',    'fn:spec');
%       cfgData.eeg  = guiPath_panel('trace', 'bottom', 'sleep_sig:eeg', ...
%           'height', 1.2);
%       guiPath(basepath, 'cfgData', cfgData);
%
%   8. Set the target. Exactly one eventTicks or stateStrip panel is curated.
%      To curate EDs from ed.mat:
%       cfgData.evt = guiPath_panel('eventTicks', 'top', 'ed', ...
%           'name', 'eventTicks');
%      cfgGui.save then says where the result is written (see REFERENCE).
%
%   9. Switch presets, reusing what is loaded. Pass the current (full) cfgData;
%      matching panels carry their data over, so only new panels load. Share a
%      ctx so shared files are not re-read.
%       [cfgData, cfgGui] = guiPath_presets('Ripples', cfgData);
%       ctx = guiPath_ctx(basepath, basename);
%       cfgData = guiPath_load(cfgData, basepath, ctx);
%      (Inside the viewer, the Preset dropdown does this with the session ctx.)
%
%   10. Reopen fast. guiPath returns the full cfgData; the live one (with
%       anything loaded through the GUI) is in hFig.UserData.cfgData.
%       [hFig, cfgData] = guiPath(basepath, 'preset', 'EDs');
%       guiPath(basepath, 'cfgData', cfgData);    % data already in hand
%
%   11. Load one more source while the viewer is open. Click Load..., pick a
%       Type, then a Source (Workspace, File, or Binary channel). No call.
%
%   12. Save. Ctrl+S, or the Save button, writes to cfgGui.save. Each save first
%       backs up any existing file (backup_file).
%
%
% REFERENCE
%
%   Panel types (the draw function follows the type; all live in guiPath_draw).
%   - trace        a 1-D signal.
%   - traces       a vertical stack of binary channels (a bin: source; kept in
%                  the file's native int16 and sliced per window, not averaged).
%   - spec         a spectrogram (adapter struct .s / .freq / .tstamps).
%   - hypnogram    read-only sleep-state strip (bout times).
%   - raster       spike raster (cell of spike-time vectors [s]).
%   - eventTicks   event marks: a tick strip in the Top, spanning lines across
%                  the signals in the Bottom (an overlay, so it takes no tile).
%                  A curation target (events mode); the Ops button overrides the
%                  Bottom look (lines <-> ticks).
%   - stateStrip   editable per-epoch label strip; a target (states mode).
%
%   Panel fields (guiPath_panel).
%   - type      one of the six above.
%   - region    'top' (full-session overview) | 'bottom' (moving window).
%   - src       the address, or an inline value.
%   - name      input identity; panels sharing a name share one loaded input, so
%               a hypnogram can appear top and bottom from a single load.
%   - label     y-axis label.
%   - height    relative panel height.
%   - clr       trace colour.
%   - ylim      [lo hi] absolute | p, a scalar clipping to the [p, 100-p]
%               percentile (0 <= p < 50; raise it when a trace looks thin) |
%               'prc', the default percentile, which is what a trace gets when
%               unset | 'full' or [] to autoscale.
%   After guiPath_load a panel also carries .data and .fs.
%
%   Address grammar (src), resolved by guiPath_src unless noted.
%   - 'ws:VAR[.path]'       base-workspace variable VAR, then a dotted path.
%   - 'bin:CH'              channel CH of <basename>.lfp (read binary).
%   - 'bin:CH>FILE'         channel CH of a given binary FILE.
%   - 'sleep_sig[:field]'   assembled sleep signals (ed_sigLoad); a field of the
%                           sSig struct if given, else the whole struct.
%   - 'FILE[:VAR.path]'     <basename>.FILE.mat, then VAR (or its sole
%                           variable) and a dotted path. E.g. 'ed',
%                           'sleep_states:ss.bouts.times'.
%   - 'fn:NAME[.field]'     a computed source (resolved by guiPath_load):
%                           'fn:spec' the spectrogram; 'fn:ripple.lfp' /
%                           'fn:ripple.filt' the ripple LFP and its filtered
%                           trace; 'fn:edLfp' the channel ED detection ran on.
%   - a value               a struct or array used in place of the string.
%
%   cfgGui fields (behaviour).
%   - name      display name.
%   - file      session file token (used to auto-detect a preset).
%   - mode      'events' | 'states' (derived from the target if omitted).
%   - win       window width [s].
%   - save      target token ('ed' | 'ripp' | 'labelsMan'), '', or a save(x)
%               handle. 'ed' / 'ripp' write the accepted mask into
%               <basename>.<token>.mat; 'labelsMan' writes state labels into
%               <basename>.sleep_labelsMan.mat.
%
%   Files. The family is one entry point plus its parts: guiPath is the tool,
%   every guiPath_* beside it is an internal of that tool. (Contrast guiTbl_*,
%   where each file is a separate tool.)
%   - guiPath.m          the viewer: figure, layout, navigation, save, Load.
%   - guiPath_presets.m  name -> [cfgData, cfgGui]; per-modality builders.
%   - guiPath_panel.m    builds one panel entry (with per-type defaults).
%   - guiPath_src.m      resolves one address to raw data.
%   - guiPath_load.m     fills data into a cfgData (skips loaded panels).
%   - guiPath_draw.m     draws one panel; the dispatch on panel type.
%   - guiPath_ctx.m      the session location + file cache.
%   - gui_loadDialog.m   the progressive Load dialog.
%   - gui_eventPanel.m   the accept / reject stepper (events mode).
%   - gui_statePanel.m   the state-assignment stepper (states mode).
%
%
% SEE ALSO
% - guiPath
% - guiPath_presets
% - guiPath_panel
% - guiPath_load
% - guiPath_src
% - guiPath_draw
% - guiPath_ctx
%
% HISTORY
% - 260705          created (declarative guiPath_presets / guiPath_load redesign).
% - 260706          rewritten: concepts + a ctx section, tighter walkthrough.
% - 260706          renamed to guiPath_doc; shared GUI package flattened to gui_*.
% - 260716          guiPath_curate renamed to guiPath; the draw functions
%                   extracted to guiPath_draw; ylim takes a scalar percentile.

end
