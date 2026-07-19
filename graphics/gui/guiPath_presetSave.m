function file = guiPath_presetSave(name, guiMap)

% Write a guiPath arrangement into a preset file.
%
% Save rewrites the VIEW half of a preset and never the DATA half:
% - an EXISTING preset is updated IN PLACE. Its recipes, its hand edits and its
%   local resolvers are kept; only the guiMap block is replaced. This is how a
%   preset's default arrangement is changed.
% - a NEW preset gets a data half that calls guiMap.base, so the base's
%   session-resolved recipes (detection channel, scaling, passband) resolve
%   again on every load and the view transfers between sessions. Freezing this
%   session's resolved values into the file would quietly show the wrong
%   channels on the next one.
%
% Either way the result is an ordinary preset file (see guiPath_doc). The two
% halves are marked by the DATA / VIEW banners, which is where Save cuts; a file
% without a VIEW banner cannot be updated in place. Any existing file is backed
% up (backup_file) before it is overwritten.
%
% EXAMPLES
% - guiPath_presetSave('ripp', guiMap)
%   replaces preset_ripp's arrangement, keeping its recipes.
% - guiPath_presetSave('rippTight', guiMap)
%   writes presets/preset_rippTight.m: guiMap's panels over guiMap.base's data.
%
% INPUTS
% - name            <char> preset token; must be a valid MATLAB name. An
%                   existing one is updated, a new one is created.
% - guiMap          <struct> .panels + .mode .win .save, and .base (which
%                   preset supplies the recipes) when NAME is new. Panel options
%                   matching the guiPath_panel default are left out; ones that
%                   cannot be written as a literal are dropped.
%
% OUTPUTS
% - file            <char> the written preset file.
%
% DEPENDENCIES
% - backup_file, guiPath_panel (for the per-type defaults).
%
% SEE ALSO
% - guiPath_preset, guiPath_panel, guiPath, guiPath_doc.
%
% HISTORY
% - 260719          created (Save preset in the guiPath control column).
% - 260719          an existing preset is updated in place (view half only);
%                   the inherit call is now only for a new preset.


%% ========================================================================
%  ARGUMENTS
%  ========================================================================
name = char(name);
assert(isvarname(name), 'guiPath_presetSave:name', ...
    'a preset name must be a valid MATLAB name (got "%s")', name);

folder = fullfile(fileparts(mfilename('fullpath')), 'presets');
if ~isfolder(folder), mkdir(folder); end
file = fullfile(folder, ['preset_', name, '.m']);


%% ========================================================================
%  BUILD THE TEXT
%  ========================================================================
% updating: the file's own data half stays, so the preset now supplies its own
% recipes and its base is itself. Creating: the base does.
if isfile(file)
    txt = spliceView(fileread(file), viewBlock(guiMap, name), name);
else
    base = '';
    if isfield(guiMap, 'base'), base = char(guiMap.base); end
    assert(~isempty(base), 'guiPath_presetSave:base', ...
        ['a new preset takes its recipes from another one, and this view ', ...
        'has none to name (it was opened from an explicit varMap)']);
    txt = newPreset(name, base, guiMap);
end


%% ========================================================================
%  WRITE
%  ========================================================================
backup_file(file);

fid = fopen(file, 'wt');
assert(fid > 0, 'guiPath_presetSave:open', 'cannot write %s', file);
fprintf(fid, '%s\n', txt);
fclose(fid);

rehash path         % so a new file is callable in this session

end


% =========================================================================
%  CODE GENERATION
% =========================================================================

function txt = spliceView(txt, viewL, name)
% swap a preset file's VIEW half for VIEWL, keeping everything else - the
% recipes, any hand edits, the local resolvers. The cut runs from the VIEW
% banner to the main function's closing 'end'.

lines = splitlines(string(txt));
iBan  = find(startsWith(strtrim(lines), '%  VIEW'), 1);
iEnd  = [];
if ~isempty(iBan)
    iEnd = iBan + find(strcmp(strtrim(lines(iBan + 1 : end)), 'end'), 1);
end
assert(~isempty(iEnd), 'guiPath_presetSave:noView', ...
    ['preset_%s.m has no "VIEW (guiMap)" banner above the function''s ', ...
    'end, so its arrangement cannot be replaced. Save under a new name, ', ...
    'or add the banner (copy it from any preset file).'], name);

head = lines(1 : iBan - 2);          % everything above the banner's %% rule
txt  = strjoin([head; string(viewL(:)); ""; lines(iEnd : end)], newline);
end


function txt = newPreset(name, base, guiMap)
% a whole new preset file: a data half that calls BASE, then the view half

L = [{sprintf('function [varMap, guiMap] = preset_%s(ctx)', name), ''}, ...
    header(name, base), ...
    banner('%%', 'DATA (varMap)'), ...
    {sprintf('varMap = preset_%s(ctx);', base), ''}, ...
    viewBlock(guiMap, base), ...
    {'', 'end', '', '% EOF'}];
txt = strjoin(L, newline);
end


function L = header(name, base)
% the help block: where the data comes from, and how to change either half

L = { ...
    '% guiPath preset saved from the GUI.', ...
    '%', ...
    sprintf(['%% The data half calls preset_%s, so its session-resolved ', ...
    'recipes (the'], base), ...
    '% detection channel, its scaling, the passband) resolve again on', ...
    '% every load and this view transfers between sessions. Press Save on', ...
    '% this preset''s name to replace the arrangement below; edit the data', ...
    '% half by hand to change what is loaded.', ...
    '%', ...
    '% INPUTS', ...
    ['% - ctx             <struct> var_ctx: basepath, basename, shared ', ...
    'file cache.'], ...
    '%', ...
    '% OUTPUTS', ...
    '% - varMap          <struct> name -> recipe.', ...
    '% - guiMap          <struct> .panels + .mode .win .save.', ...
    '%', ...
    '% SEE ALSO', ...
    '% - guiPath_preset, guiPath_presetSave, guiPath_panel, guiPath_doc.', ...
    '%', ...
    '% HISTORY', ...
    sprintf('%% - %s          saved from guiPath as "%s".', stamp(), name), ...
    '', ''};
end


function L = viewBlock(guiMap, base)
% the VIEW half, banner included: the behaviour struct then one line per panel,
% in stacking order

L = banner('%%', 'VIEW (guiMap)');

opts = [{sprintf('%s, %s', lit('base'), lit(base))}, ...
    behaviourOpt(guiMap, 'mode'), behaviourOpt(guiMap, 'win'), ...
    behaviourOpt(guiMap, 'save')];
L{end + 1} = wrapCall('guiMap = struct(''panels'', struct()', opts);

fns = fieldnames(guiMap.panels);
for iPan = 1 : numel(fns)
    L{end + 1} = panelLine(fns{iPan}, guiMap.panels.(fns{iPan})); %#ok<AGROW>
end
end


function L = banner(lead, title)
% a foldable section banner, as three lines
rule = repmat('=', 1, 72);
L = {[lead, ' ', rule], ['%  ', title], ['%  ', rule]};
end


function opt = behaviourOpt(guiMap, fld)
% one 'name, value' pair for the guiMap struct call, or nothing when the field
% is absent or not writable as a literal (a save handle, say)

opt = {};
if ~isfield(guiMap, fld) || isempty(guiMap.(fld)), return; end
s = lit(guiMap.(fld));
if isempty(s), return; end
opt = {sprintf('%s, %s', lit(fld), s)};
end


function line = panelLine(field, pan)
% one 'guiMap.panels.<field> = guiPath_panel(...)' statement. Only options that
% differ from the type's default are written, so the file stays readable.

head = sprintf('guiMap.panels.%s = guiPath_panel(%s, %s, %s', field, ...
    lit(pan.type), lit(pan.region), lit(pan.var));
dflt = guiPath_panel(pan.type, pan.region, pan.var);

opts = {};
optFns = setdiff(fieldnames(dflt), {'type', 'region', 'var'}, 'stable');
for iFld = 1 : numel(optFns)
    fld = optFns{iFld};
    if ~isfield(pan, fld) || isequaln(pan.(fld), dflt.(fld)), continue; end
    s = lit(pan.(fld));
    if isempty(s), continue; end        % not a literal -> keep the default
    opts{end + 1} = sprintf('%s, %s', lit(fld), s); %#ok<AGROW>
end

line = wrapCall(head, opts);
end


function line = wrapCall(head, opts)
% close a call, keeping every line under 80 columns: options go on continuation
% lines indented by 4, packed as tightly as they fit

if isempty(opts)
    line = [head, ');'];
    return
end

line = [head, ', ', strjoin(opts, ', '), ');'];
if numel(line) <= 80, return; end

line = head;
cur  = '';
for iOpt = 1 : numel(opts)
    if isempty(cur), cand = ['    ', opts{iOpt}];
    else,            cand = [cur, ', ', opts{iOpt}];
    end
    if numel(cand) > 76 && ~isempty(cur)
        line = sprintf('%s, ...\n%s', line, cur);
        cur  = ['    ', opts{iOpt}];
    else
        cur = cand;
    end
end
line = sprintf('%s, ...\n%s);', line, cur);
end


function s = lit(v)
% a MATLAB literal for a panel / behaviour value; '' when it cannot be written
% (a function handle, a cell, an oversized array)

s = '';
if ischar(v)
    s = ['''', strrep(v, '''', ''''''), ''''];
elseif isstring(v) && isscalar(v)
    s = lit(char(v));
elseif islogical(v) || (isnumeric(v) && isreal(v) && numel(v) <= 8)
    s = mat2str(v);
end
end


function s = stamp()
% today as YYMMDD, for the HISTORY line
s = char(datetime('now', 'Format', 'yyMMdd'));
end

% EOF
