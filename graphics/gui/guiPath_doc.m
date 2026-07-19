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
% Everything rests on one split: the DATA and the VIEW are two separate maps.
%
%
% CONCEPTS
%
%   A preset is two maps, held and passed separately.
%   - varMap    the data: a struct, one field per signal / set, each field a
%               recipe (see var_recipe) saying WHERE the data lives. io/var_load
%               fills each recipe with its .data / .fs, in place. This is the
%               reusable I/O layer - the same var_load / var_fetch serve any
%               pipeline, not just the GUI.
%   - guiMap    the arrangement: .panels (one view panel per field, see
%               guiPath_panel) plus behaviour (.name .base .mode .win .save). A
%               panel says WHAT it is (type), WHERE it sits (region) and WHICH
%               data it draws (var, a varMap field). This is what you save.
%   guiPath_preset(name, basepath) returns both. guiPath draws them.
%
%   Data and view are joined by name. A panel's .var picks a varMap field;
%   SEVERAL panels may name one field. That is how a hypnogram shows top and
%   bottom from a single load, and how an event set shows as Top ticks plus a
%   Bottom overlay - one datum, many views, no duplication.
%
%   A recipe is an indirection. A panel does not hold its data or know how to
%   read it; the varMap field it names holds a recipe - a small struct with a
%   kind (matvar | matfield | bin | ws | value) and, optionally, a transform
%   chain and a dot-path. var_fetch is the single place that turns a recipe into
%   raw data (see the grammar in REFERENCE). Session-specific choices (which
%   channel ripple / ED detection ran on) are resolved by the preset and frozen
%   into plain recipes, so the loader stays generic.
%
%   Slim versus full. A varMap from guiPath_presets is slim: recipes only, no
%   data. var_load fills each field's .data / .fs, returning the same struct now
%   full. A field that already carries data is skipped, so topping up is cheap.
%
%   Loader is view-blind; the view shapes. var_load / var_fetch return the raw
%   value. guiPath_shape then turns it into what a panel draws: an events struct,
%   an hours-cell hypnogram, a channel stack's stats, a resolved y-limit. This
%   happens once, when a preset lands (see mapsToConfig in guiPath).
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
% CONCEPTS: a preset is a file
%
%   Each preset is one MATLAB file in graphics/gui/presets, preset_<token>.m,
%   that returns the two maps for a session:
%
%       function [varMap, guiMap] = preset_ripp(ctx)
%
%   The file IS the preset. Its token is what the PRESET selector shows, what
%   'preset' matches, and - when <basename>.<token>.mat exists - what
%   auto-detection picks. Adding a preset means adding a file; nothing else
%   knows the list (guiPath_preset scans the folder).
%
%   Every preset file has two halves, marked by banners:
%       DATA (varMap)   var_recipe calls - WHAT is loaded.
%       VIEW (guiMap)   guiPath_panel calls - HOW it is shown.
%   A guiMap always needs a varMap: a panel names a var, and a panel naming a
%   var the varMap does not have is dropped when the view is built. That is why
%   a preset returns BOTH halves, and why an arrangement cannot stand alone.
%
%   Save rewrites the VIEW half, never the DATA half. Arrange a session in the
%   GUI (add / delete / reorder panels, set an amplitude, switch a Bottom event
%   set to ticks), press Save next to the PRESET selector, and name it:
%   - a name already in the list UPDATES that preset in place. Its recipes,
%     hand edits and local resolvers are kept; only the guiMap block is
%     replaced. This is how a modality's default arrangement is changed - press
%     Save, type ripp.
%   - a new name CREATES a preset whose data half calls the current one:
%
%       varMap = preset_ripp(ctx);
%       guiMap = struct('panels', struct(), 'base', 'ripp', ...);
%       guiMap.panels.spec = guiPath_panel('spec', 'top', 'spec');
%
%   Why a new preset calls rather than copies. Ripple and ED presets resolve
%   session-specific values - the detection channel, its bit2uv scaling, the
%   passband - when they are built. Writing THIS session's resolved numbers into
%   the file would quietly show the wrong channels on the next session, so the
%   file re-runs its base instead. Same reason for a trace's y-limits: what is
%   saved is the declared 'prc' / 'full', not the numbers this session resolved.
%
%   Two things to know.
%   - To change what a preset LOADS, edit the DATA half by hand; Save will not
%     touch it. In a created preset that means replacing the one call with your
%     own var_recipe lines.
%   - A panel over something added with Load... has no recipe in the preset, so
%     it cannot be written; Save says which ones will be missing.
%
%   Renaming a preset means renaming its FILE and its function line together.
%   guiPath_preset refuses a file whose declared name does not match it - a
%   half-rename would otherwise make a created preset call itself.
%
%
% CONCEPTS: ctx, the session cache
%
%   ctx is a small struct that travels with a load. var_ctx(basepath, basename)
%   builds it. It carries the session location (.basepath, .basename) and one
%   cache (.cache) shared by every recipe resolved during the load.
%
%   The problem it solves. Many fields read from the same file. In the EDs preset
%   the hypnogram, spectrogram, EMG and EMG RMS all come out of the assembled
%   sleep signals (sleep_sig.mat, an expensive read). The cache stores each file
%   the first time it is read and returns the stored copy on every later hit, so
%   a given file is read at most once per ctx.
%
%   The mechanism (worth understanding, because MATLAB makes it surprising). A
%   struct is a value: passing it into a function copies it, and writes inside
%   never reach the caller. containers.Map is one of MATLAB's few reference
%   (handle) types: copying the ctx struct copies the handle, not the map, so
%   every copy points at the same underlying map. A write through any copy
%   (ctx.cache(key) = value) is seen through all of them. That is why var_fetch
%   can take ctx by value, write into ctx.cache, never return ctx, and still
%   leave the caller's cache filled.
%
%   Two levels of "load once" stack. Field level: var_load skips a field that
%   already has data. File level: even when a field must load, ctx.cache makes
%   the underlying file read happen once, so two fields into one file
%   (sleep_sig emg and emg_rms) share the one cached read.
%
%   Use case, a preset switch. guiPath builds one ctx when it opens and keeps it
%   for the session (hFig.UserData.ctx). Switching EDs -> Ripples builds the new
%   preset's varMap and var_loads it with the SAME ctx, so the shared files
%   (sleep_sig, session, spikes) are already cached and only the ripple-specific
%   reads run.
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
%       guiPath(basepath, 'preset', 'sleep_states');  % force one by token
%
%   2. Get a preset's data + arrangement.
%       [varMap, guiMap] = guiPath_preset('ed', basepath);   % ripp sleep_states
%       fieldnames(varMap)'          % the signals / sets
%       fieldnames(guiMap.panels)'   % the panels, in stacking order
%       guiPath_preset()             % the available tokens
%       edit preset_ed               % the file itself
%
%   3. Read one field's recipe and one panel.
%       varMap.spec                  % a matfield recipe (the spectrogram)
%       guiMap.panels.spec.type      % 'spec'
%       guiMap.panels.spec.region    % 'top'
%       guiMap.panels.spec.var       % 'spec' (the varMap field it draws)
%
%   4. Load the data. Every recipe is filled with .data; a field that already
%      has data is left alone.
%       varMap = var_load(varMap, basepath);
%       varMap.spec.data             % now present
%
%   5. Open your own maps. A slim varMap (no data) is loaded for you; omit guiMap
%      and it is derived.
%       guiPath(basepath, 'varMap', varMap, 'guiMap', guiMap);
%
%   6. Edit or add a panel + its data. guiPath_panel(type, region, var) builds a
%      view panel; its field position sets its stacking order.
%       guiMap.panels.emg.height = 1.0;
%       varMap.emg2 = var_recipe('matfield', 'file', 'sleep_sig', 'field', 'emg');
%       guiMap.panels.emg2 = guiPath_panel('trace', 'bottom', 'emg2', 'label', 'EMG 2');
%       guiPath(basepath, 'varMap', varMap, 'guiMap', guiMap);
%
%   7. Build from a blank template.
%       varMap = struct('spec', var_recipe('matfield', 'file', 'sleep_sig', ...
%           'field', {{'spec', 'spec_freq', 'spec_tstamps'}}), ...
%           'eeg', var_recipe('matfield', 'file', 'sleep_sig', 'field', 'eeg'));
%       guiMap = struct('panels', struct( ...
%           'spec', guiPath_panel('spec',  'top',    'spec'), ...
%           'eeg',  guiPath_panel('trace', 'bottom', 'eeg', 'height', 1.2)));
%       guiPath(basepath, 'varMap', varMap, 'guiMap', guiMap);
%
%   8. Set the target. Exactly one eventTicks or stateStrip panel is curated.
%      To curate EDs from ed.mat:
%       varMap.ed = var_recipe('matvar', 'file', 'ed');
%       guiMap.panels.evt = guiPath_panel('eventTicks', 'top', 'ed');
%      guiMap.save then says where the result is written (see REFERENCE).
%
%   9. Reopen fast. guiPath returns the loaded varMap; the live one (with
%      anything loaded through the GUI) is in hFig.UserData.varMap / .guiMap.
%       [hFig, varMap, guiMap] = guiPath(basepath, 'preset', 'ed');
%       guiPath(basepath, 'varMap', varMap, 'guiMap', guiMap);  % data in hand
%
%   10. Load one more source while the viewer is open. Click Load..., pick a
%       Type, then a Source (Workspace, File, or Binary channel). No call.
%
%   11. Keep an arrangement. Arrange the panels, press Save beside the PRESET
%       selector, and name it: an existing name replaces that preset's
%       arrangement, a new name creates one over the same data (see CONCEPTS: a
%       preset is a file). Either is then in the selector, for every session.
%       guiPath_presetSave('ripp', guiMap);        % the same, from code
%
%   12. Save curation. Ctrl+S, or the action panel's Save, writes to
%       guiMap.save. Each save first backs up any existing file (backup_file).
%
%
% REFERENCE
%
%   Panel types (the draw function follows the type; all live in guiPath_draw).
%   - trace        a 1-D signal.
%   - traces       a vertical stack of binary channels (a bin recipe kept native
%                  int16 and sliced per window, not averaged).
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
%   - type      one of the seven above.
%   - region    'top' (full-session overview) | 'bottom' (moving window).
%   - var       the varMap field this panel draws. Panels sharing a var share
%               one loaded input, so a hypnogram can appear top and bottom.
%   - label     y-axis label.
%   - height    relative panel height.
%   - clr       trace / tick colour.
%   - ylim      [lo hi] absolute | p, a scalar clipping to the [p, 100-p]
%               percentile (0 <= p < 50; raise it when a trace looks thin) |
%               'prc', the default percentile, which a trace gets when unset |
%               'full' or [] to autoscale.
%   - render    how a type that has a choice draws HERE: a Bottom event set is
%               'lines' (spanning, no tile) or 'strip' (a tick lane). '' = by
%               type + region. This is what the Ops button sets.
%   - yAdjust   amplitude factor for a trace / traces / spec, as shift+scroll
%               sets it live. 1 = as loaded.
%
%   Recipe grammar (var_recipe). Resolution is fetch (by kind) -> transform
%   chain -> dot-path.
%   - matvar    a field of <basename>.file.mat via its wrapper variable, then a
%               path. var_recipe('matvar','file','sleep_states','var','ss', ...
%               'path','bouts.times').
%   - matfield  a named top-level field of a -struct .mat (e.g. sleep_sig eeg /
%               emg / emg_rms / spec). A cellstr field packs several into a
%               struct. var_recipe('matfield','file','sleep_sig','field','eeg').
%   - bin       binary channel(s); nCh / fs from session.mat. 'average' means to
%               a single trace; 'outClass','native' keeps int16 for a stack.
%               var_recipe('bin','file','lfp','ch',[5 6 7],'average',false, ...
%               'outClass','native').
%   - ws        a base-workspace variable + path (the Load dialog).
%   - value     an already-materialized value, used inline.
%   - transform ops run after the fetch: eegSub, emg, emgRms, spec, rippPrep.
%               var_recipe('bin','file','lfp','ch',ch, ...
%               'transform',{{'rippPrep',{[80 250]}}},'path','filt').
%   A bare 'file.path' string is shorthand for a matvar read.
%
%   guiMap behaviour.
%   - name      the preset token; set by guiPath_preset from the file name.
%   - base      which preset supplies the recipes: itself for a preset whose
%               data half is its own, the called token for a created one, ''
%               for a hand-built varMap (which can update an existing preset's
%               arrangement, but cannot create a new preset).
%   - mode      'events' | 'states' (derived from the target if omitted).
%   - win       window width [s].
%   - save      target token ('ed' | 'ripp' | 'labelsMan'), '', or a save(x)
%               handle. 'ed' / 'ripp' write the accepted mask into
%               <basename>.<token>.mat; 'labelsMan' writes state labels into
%               <basename>.sleep_labelsMan.mat.
%
%   Files. The family is one entry point plus its parts. The data layer (io/) is
%   shared, not GUI-specific.
%   - guiPath.m          the viewer: figure, layout, navigation, save, Load.
%   - presets/           one preset_<token>.m per preset; the whole list.
%   - guiPath_preset.m   finds and calls one; no arg -> the available tokens.
%   - guiPath_presetSave.m  writes a live arrangement into a preset file.
%   - guiPath_panel.m    builds one view panel (with per-type defaults).
%   - guiPath_shape.m    shapes a raw value into a drawn input, by type.
%   - guiPath_draw.m     draws one panel; the dispatch on panel type.
%   - io/var_recipe.m    builds one data recipe.
%   - io/var_fetch.m     resolves one recipe to its raw value.
%   - io/var_load.m      fills a varMap of recipes, in place.
%   - io/var_ctx.m       the session location + file cache.
%   - gui_loadDialog.m   the progressive Load dialog.
%   - gui_eventPanel.m   the accept / reject stepper (events mode).
%   - gui_statePanel.m   the state-assignment stepper (states mode).
%
%
% SEE ALSO
% - guiPath
% - guiPath_preset
% - guiPath_presetSave
% - guiPath_panel
% - guiPath_shape
% - var_recipe
% - var_load
% - guiPath_draw
%
% HISTORY
% - 260705          created (declarative guiPath_presets / guiPath_load redesign).
% - 260706          rewritten: concepts + a ctx section, tighter walkthrough.
% - 260716          guiPath_curate renamed to guiPath; the draw functions
%                   extracted to guiPath_draw; ylim takes a scalar percentile.
% - 260719          data / view split onto one reusable loader: varMap (recipes,
%                   io/var_load) + guiMap (view panels); the address grammar and
%                   guiPath_load / guiPath_src / guiPath_ctx are gone.
% - 260719          presets became one file each (presets/preset_<token>.m,
%                   found by guiPath_preset); Save writes the live arrangement
%                   over a base preset's data (guiPath_presetSave).

end
