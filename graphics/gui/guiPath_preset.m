function [varMap, guiMap] = guiPath_preset(name, ctx)

% Find and load a guiPath preset.
%
% A preset is ONE MATLAB file in graphics/gui/presets, named preset_<name>.m,
% that returns a session's two halves:
%
%   function [varMap, guiMap] = preset_ripp(ctx)
%
% The file IS the preset. Its token (<name>) is the identity shown in the GUI
% and, when <basename>.<name>.mat exists, what auto-detection matches. Each file
% has a DATA half (var_recipe calls) and a VIEW half (guiPath_panel calls); edit
% the first by hand, and guiPath_presetSave rewrites the second from a view
% arranged in the GUI.
%
% Renaming a preset means renaming its file AND its function line. A created
% preset takes its recipes by calling another preset by name, so a file renamed
% to the one it calls would call itself; this function refuses to run a file
% whose declared name does not match it.
%
% This function only finds and calls the file, then stamps .name (the token)
% and .base (which preset supplies the recipes - itself, unless the file says
% otherwise). Missing behaviour fields are filled by guiPath when it opens.
%
% EXAMPLES
% - [varMap, guiMap] = guiPath_preset('ripp', basepath)
%   the ripple preset for a session; pass varMap to var_load to fill it.
% - names = guiPath_preset()
%   every available token, for the dropdown and auto-detection.
%
% INPUTS
% - name            <char>(opt) preset token. No arg -> the available tokens.
% - ctx             <struct | char>(opt) var_ctx, whose file cache is shared
%                   with the later var_load, or a basepath. Default pwd.
%
% OUTPUTS
% - varMap          <struct> name -> recipe (see var_recipe).
% - guiMap          <struct> .panels + .name .base (+ what the preset set).
% - names           <cellstr> the available tokens (the no-arg form).
%
% DEPENDENCIES
% - var_ctx, and the preset files in ./presets.
%
% SEE ALSO
% - guiPath_presetSave, var_recipe, guiPath_panel, guiPath, guiPath_doc.
%
% HISTORY
% - 260719          created; presets became one file each (were builders inside
%                   guiPath_presets).


%% ========================================================================
%  ENUMERATION (no arg -> the available tokens)
%  ========================================================================
folder = fullfile(fileparts(mfilename('fullpath')), 'presets');

if nargin < 1 || isempty(name)
    dirFiles = dir(fullfile(folder, 'preset_*.m'));
    varMap = cell(1, numel(dirFiles));
    for iFile = 1 : numel(dirFiles)
        [~, stem] = fileparts(dirFiles(iFile).name);
        varMap{iFile} = stem(numel('preset_') + 1 : end);
    end
    varMap = sort(varMap);
    return
end


%% ========================================================================
%  LOAD ONE
%  ========================================================================
if nargin < 2 || isempty(ctx), ctx = pwd; end
if ~isstruct(ctx), ctx = var_ctx(char(ctx)); end

name = char(name);
fcn  = ['preset_', name];
file = fullfile(folder, [fcn, '.m']);
if ~isfile(file)
    error('guiPath_preset:name', 'no preset "%s" in %s', name, folder);
end

% a preset is identified by its FILE name, so the function it declares must
% match. Renaming the file alone leaves MATLAB calling it by the file name while
% its body still calls the old name - and since a saved preset inherits by
% calling its base, a renamed one calls ITSELF: unbounded recursion, out of
% memory, no useful error. Catch it here instead.
declared = declaredFcn(file);
if ~isempty(declared) && ~strcmp(declared, fcn)
    error('guiPath_preset:renamed', ['%s.m declares "%s". A preset is ', ...
        'named by its file, so rename the function line to "%s" too (or ', ...
        'rename the file back to %s.m).'], fcn, declared, fcn, declared);
end

if isempty(which(fcn))
    addpath(folder);            % folder not on the path (or a just-saved file)
end

[varMap, guiMap] = feval(fcn, ctx);

if ~isstruct(guiMap), guiMap = struct(); end
guiMap.name = name;
if ~isfield(guiMap, 'base') || isempty(guiMap.base)
    guiMap.base = name;         % a preset supplies its own recipes
end

end


% =========================================================================
%  HELPERS
% =========================================================================

function nm = declaredFcn(file)
% the function name a .m file declares ('' when it declares none). Anchored to a
% line start, so a 'function ...' line inside the help block never matches.

nm = '';
tok = regexp(fileread(file), '^\s*function\s+(?:[^=\n]*=\s*)?([A-Za-z]\w*)', ...
    'tokens', 'once', 'lineanchors');
if ~isempty(tok), nm = tok{1}; end
end

% EOF
