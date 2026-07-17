function varMap = varMap_build(spec, varargin)

% Build or complete a varMap's view, and optionally load its data.
%
% The viewer's half of the layer: it owns the presets and the view. var_load
% owns the data and never touches view. varMap_build may call var_load;
% var_load never calls varMap_build.
%
% A view is added in three tiers. A preset writes it by hand. Otherwise the
% type is taken from the recipe (fn -> its panel, args -> a trace) and the rest
% from a per-type table. A plain .mat read cannot be typed from the recipe
% alone, so it is left view-less rather than guessed.
%
% EXAMPLES
% - varMap = varMap_build('ed', 'basepath',bp, 'flgLoad',true)
%   builds the 'ed' preset (recipes and their view) and loads the data.
% - varMap = varMap_build(myMap)
%   takes a partial map (recipes only) and adds a view to each entry.
%
% INPUTS
% - spec            <char | struct> a preset name, or a partial varMap.
%
% NAME-VALUE
% - basepath        <char> session folder. Default pwd.
% - flgLoad         <logical> also run var_load at the end. Default false.
%
% OUTPUTS
% - varMap          <struct> recipes + view (+ data if flgLoad).
%
% DEPENDENCIES
% - var_load (only when flgLoad is true).
%
% HISTORY
% - 260717          created (unified var_* I/O layer).


%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addParameter(p, 'basepath', pwd);
addParameter(p, 'flgLoad',  false, @islogical);
parse(p, varargin{:});

basepath = p.Results.basepath;
flgLoad  = p.Results.flgLoad;


%% ========================================================================
%  BUILD
%  ========================================================================

if ischar(spec) || isstring(spec)
    varMap = preset(char(spec));            % a named layout
else
    varMap = spec;                          % a user's partial map
end

varMap = addView(varMap);                   % fill any missing view

if flgLoad
    varMap = var_load(varMap, basepath);
end

end


% =========================================================================
%  HELPERS
% =========================================================================

function varMap = preset(name)
% named layouts. Channels are placeholders; a real preset reads them from
% session.channelTags. Add a preset -> add a case.
switch name
    case 'ed'
        varMap = struct();
        varMap.spec = struct('file', 'lfp', 'args', {{'ch', 1}}, ...
            'fn', 'spec', 'view', viewDefaults('spec'));
        varMap.eeg  = struct('file', 'lfp', 'args', {{'ch', 1}}, ...
            'view', viewDefaults('trace'));
        varMap.ed   = struct('file', 'ed', ...
            'view', viewDefaults('eventTicks'));

    otherwise
        error('varMap_build:preset', 'unknown preset "%s"', name);
end
end


function varMap = addView(varMap)
% give every view-less entry a default view, when its type can be told
fldNames = fieldnames(varMap);
for iFld = 1 : numel(fldNames)
    entry = varMap.(fldNames{iFld});
    if isstruct(entry) && isfield(entry, 'view')
        continue                            % view already set
    end
    type = inferType(entry);
    if isempty(type)
        continue                            % cannot type it -> leave view-less
    end
    entry.view = viewDefaults(type);
    varMap.(fldNames{iFld}) = entry;
end
end


function type = inferType(entry)
% the panel type a recipe implies, or '' when it cannot be told
if ~isstruct(entry)
    type = '';                              % a bare-string read
elseif isfield(entry, 'fn')
    switch entry.fn
        case 'spec',   type = 'spec';
        case 'emgRms', type = 'trace';
        otherwise,     type = '';
    end
elseif isfield(entry, 'args')
    type = 'trace';                         % a binary channel
else
    type = '';                              % a .mat read
end
end


function view = viewDefaults(type)
% per-type layout: where it sits, how tall, its stacking order
switch type
    case 'spec'
        view = struct('type', 'spec', 'region', 'top', ...
            'height', 1.4, 'order', 20, 'label', 'Freq (Hz)');
    case 'trace'
        view = struct('type', 'trace', 'region', 'bottom', ...
            'height', 1.0, 'order', 50, 'label', '');
    case 'traces'
        view = struct('type', 'traces', 'region', 'bottom', ...
            'height', 2.5, 'order', 50, 'label', 'LFP');
    case 'eventTicks'
        view = struct('type', 'eventTicks', 'region', 'top', ...
            'height', 0.28, 'order', 40, 'label', 'Events');
    case 'stateStrip'
        view = struct('type', 'stateStrip', 'region', 'top', ...
            'height', 0.5, 'order', 15, 'label', 'State');
    case 'hypnogram'
        view = struct('type', 'hypnogram', 'region', 'top', ...
            'height', 0.28, 'order', 10, 'label', 'State');
    case 'raster'
        view = struct('type', 'raster', 'region', 'bottom', ...
            'height', 1.2, 'order', 60, 'label', 'Units');
    otherwise
        view = struct('type', type, 'region', 'bottom', ...
            'height', 1.0, 'order', 99, 'label', '');
end
end

% EOF
