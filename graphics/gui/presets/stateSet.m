function sSet = stateSet(ctx)

% Composes the sleep-state set every preset shows and curates.
%
% One per-epoch label vector on the session's own state config, ready for a
% guiPath 'stateStrip' panel. It sits with the presets rather than in the
% guiPath_* framework because it is the one place that knows the AccuSleep
% layout; the framework stays modality-agnostic.
%
% The labels are ss.labels - the classifier's output with manual scoring
% already merged in (as_classify) - with <basename>.sleep_labelsMan.mat merged
% on top again, so scoring saved from the GUI shows on the next open, before
% as_classify has re-run. Reading labelsMan alone shows only the epochs a human
% touched, which on a classified session is a nearly empty strip.
%
% Names and colours come from ss.info, so the strip matches the states this
% session was actually scored with. The assignable range is the NAME count, not
% cfg.nstates: as_loadConfig counts 6 states but ships 7 names, BIN being the
% 7th, and undefined is one past the last name (AccuSleep's unscored code).
%
% INPUTS
% - ctx             <struct> var_ctx: basepath, basename, shared file cache.
%
% OUTPUTS
% - sSet            <struct> .labels .epochT .names .colors, or [] when the
%                   session holds no sleep scoring.
%
% SEE ALSO
% - guiPath_panel, guiPath_draw, var_recipe, as_classify.
%
% HISTORY
% - 260720          created; one state set for every preset (ripp / ed drew a
%                   read-only hypnogram from ss.bouts.times while the states
%                   preset composed its own strip from sleep_labelsMan).
% - 260720b         returns [] instead of throwing on an unscored session. A
%                   preset resolves its session values while it is BUILT, not
%                   while it loads, so an error here took the whole preset down
%                   before var_load could drop the entry.

% A session that was never sleep-scored has no state set to show. Returning
% empty rather than throwing is what lets a preset open on it at all: the
% caller omits the varMap entry and guiPath drops the panels that name it, so
% an unscored recording (the EA cohort has no sleep_states) still opens on its
% signals and its events.
sSet = [];
try
    ss = var_fetch(var_recipe('matvar', 'file', 'sleep_states', ...
        'var', 'ss'), ctx);
    epochT = var_fetch(var_recipe('matfield', 'file', 'sleep_sig', ...
        'field', 'spec_tstamps'), ctx);
catch
    return
end
if ~isstruct(ss) || ~isfield(ss, 'info') || ~isfield(ss, 'labels')
    return
end

names  = ss.info.names(:)';
colors = ss.info.colors(:)';
labels = double(ss.labels(:));

% the manual layer overrides wherever it is scored - the same merge
% as_classify makes when it rebuilds ss
man = manLabels(ctx, numel(labels));
if ~isempty(man)
    iMan = man <= numel(names);
    labels(iMan) = man(iMan);
end

sSet = struct('labels', labels, 'epochT', epochT(:), ...
    'names', {names}, 'colors', {colors});

end


% =========================================================================
%  SESSION RESOLVERS
% =========================================================================

function man = manLabels(ctx, nEp)
% <basename>.sleep_labelsMan.mat, or [] when it is absent or its length has
% drifted from the epoch count (a positional vector cannot be remapped)

man = [];
file = fullfile(ctx.basepath, [ctx.basename, '.sleep_labelsMan.mat']);
if ~isfile(file), return; end
s = load(file, 'labels');
if isfield(s, 'labels') && numel(s.labels) == nEp
    man = double(s.labels(:));
end
end

% EOF
