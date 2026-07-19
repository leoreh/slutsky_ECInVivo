function p = guiPath_panel(type, region, var, varargin)

% Build one view panel for a guiPath guiMap.
%
% A guiMap is the arrangement half of a preset: a flat struct with one panel per
% field (field order = stacking order within a region). A panel says WHAT it is
% (type), WHERE it sits (region) and WHICH data it draws (var, the name of a
% varMap entry) - never the data itself. Several panels may name one var (a
% top + bottom pair shares a single loaded input); that is how one signal shows
% in two places without loading it twice.
%
% The data half is a varMap (see var_recipe); var_load fills it. The two are
% joined by name at open time (see guiPath > mapsToConfig).
%
% INPUTS
% - type            <char> 'trace' | 'traces' | 'spec' | 'hypnogram' |
%                   'raster' | 'eventTicks' | 'stateStrip'.
% - region          <char> 'top' | 'bottom' (aliases 'wide' | 'narrow').
% - var             <char> the varMap entry this panel draws.
%
% NAME-VALUE (optional; override the per-type default)
% - height          <num>  relative panel height.
% - label           <char> panel y-label.
% - clr             <ColorSpec> trace / tick colour.
% - ylim            <2-vec | scalar | 'prc' | 'full'> y-limits. A 2-vec is
%                   absolute; a scalar p clips to the [p, 100-p] percentile
%                   (0 <= p < 50); 'prc' is the default percentile (what a trace
%                   gets unset); 'full' / [] autoscale. See guiPath_shape.
% - render          <char> how the type draws here, when it has a choice: a
%                   Bottom event set is 'lines' (spanning, no tile of its own)
%                   or 'strip' (a tick lane). '' (default) = by type + region.
% - yAdjust         <num>  amplitude factor for trace / traces / spec, as set
%                   live by shift+scroll. Default 1.
%
% OUTPUTS
% - p               <struct> .type .region .var .height .label .clr .ylim
%                   .render .yAdjust.
%
% SEE ALSO
% - var_recipe, guiPath_preset, guiPath_shape, guiPath, guiPath_doc.
%
% HISTORY
% - 260719          view panel over a varMap entry (was a fused src+view+data
%                   panel; recipe moved to var_recipe, loading to var_load).
% - 260719          render + yAdjust added, so a saved preset keeps a Bottom
%                   set's lines / ticks choice and the amplitude it was given.

td = typeDefaults(type);
p = struct('type', type, 'region', region, 'var', var, ...
    'height', td.height, 'label', td.label, 'clr', 'k', 'ylim', [], ...
    'render', '', 'yAdjust', 1);
if strcmp(type, 'trace'), p.ylim = 'prc'; end     % traces auto-percentile unless set

for iArg = 1 : 2 : numel(varargin)
    key = varargin{iArg};
    v   = varargin{iArg + 1};
    switch lower(key)
        case 'height',  p.height  = v;
        case 'label',   p.label   = v;
        case 'clr',     p.clr     = v;
        case 'ylim',    p.ylim    = v;
        case 'render',  p.render  = v;
        case 'yadjust', p.yAdjust = v;
        otherwise, error('guiPath_panel:arg', 'unknown option "%s"', key);
    end
end

end


function td = typeDefaults(type)
% per-panel-type fallback height / label
switch type
    case 'spec',       td = struct('height', 1.4,  'label', 'Freq (Hz)');
    case 'hypnogram',  td = struct('height', 0.28, 'label', 'State');
    case 'eventTicks', td = struct('height', 0.28, 'label', 'Events');
    case 'stateStrip', td = struct('height', 0.5,  'label', 'State');
    case 'raster',     td = struct('height', 1.2,  'label', 'Units');
    case 'trace',      td = struct('height', 1.0,  'label', '');
    case 'traces',     td = struct('height', 2.5,  'label', 'LFP');
    otherwise,         td = struct('height', 1.0,  'label', '');
end
end

% EOF
