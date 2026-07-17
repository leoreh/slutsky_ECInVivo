function p = guiPath_panel(type, region, src, varargin)
% GUIPATH_PANEL Build one panel entry for a guiPath cfgData.
%
%   p = GUIPATH_PANEL(type, region, src, Name, Value, ...) returns a struct that
%   declares a single panel: WHAT it is (type), WHERE it sits (region) and WHERE
%   its data lives (src, an address for guiPath_src). A cfgData is a flat struct
%   with one such entry per field (field name = panel name); the field order is
%   the stacking order within a region. The per-type appearance defaults (height,
%   label, order) are filled here, so a panel is self-describing; only the data
%   is added later, by guiPath_load.
%
%   INPUTS:
%       type    - (Char) 'trace' | 'traces' | 'spec' | 'hypnogram' | 'raster' |
%                        'eventTicks' | 'stateStrip'. trace is one signal;
%                        traces is a vertical stack of binary channels (src must
%                        be a bin: address). eventTicks/stateStrip mark the
%                        curation target (they set the mode; there is at most one).
%       region  - (Char) 'top' | 'bottom' (aliases: 'wide' | 'narrow').
%       src     - (Char | value) address for guiPath_src, a computed 'fn:...'
%                 source, or an inline value.
%       Name/Value (all optional; override the per-type default):
%           'name'   (Char)   input identity; panels that share a name share one
%                             loaded input. Defaults to the cfgData field name.
%           'fs'     (Num)    sampling rate [Hz] for a trace whose address does
%                             not carry one (ws / bin sources).
%           'height' (Num)    relative panel height.
%           'clr'    (ColorSpec) trace colour.
%           'ylim'   (2-vec | scalar | 'prc' | 'full') y-limits. A 2-vec is
%                             absolute. A scalar p clips to the [p, 100-p]
%                             percentile (0 <= p < 50); raise it when a trace
%                             looks thin. 'prc' takes the default percentile
%                             and is what traces get when unset. 'full' / []
%                             autoscale. See resolveYlim in guiPath_load.
%           'label'  (Char)   panel y-label.
%           'order'  (Num)    override stacking order within the region.
%
%   OUTPUT:
%       p       - (Struct) one panel entry (guiPath_load fills .data / .fs later).
%
%   See also guiPath_load, guiPath_src, guiPath_presets, guiPath, guiPath_doc.
%
%   HISTORY:
%       Created: 05 Jul 2026 - declarative redesign (panel constructor).
%       Updated: 05 Jul 2026 - per-type defaults folded in (was curate_typeDefaults).

td = typeDefaults(type);
p = struct('type', type, 'region', region, 'src', {src}, 'name', '', ...
    'fs', [], 'height', td.height, 'clr', 'k', 'ylim', [], 'label', td.label, 'order', td.order);
if strcmp(type, 'trace'), p.ylim = 'prc'; end     % traces auto-percentile unless overridden

for i = 1:2:numel(varargin)
    key = varargin{i}; v = varargin{i + 1};
    switch lower(key)
        case 'name',   p.name   = v;
        case 'fs',     p.fs     = v;
        case 'height', p.height = v;
        case 'clr',    p.clr    = v;
        case 'ylim',   p.ylim   = v;
        case 'label',  p.label  = v;
        case 'order',  p.order  = v;
        otherwise, error('guiPath_panel:arg', 'unknown option "%s"', key);
    end
end
end

function td = typeDefaults(type)
% per-panel-type fallback height / label / stacking order
switch type
    case 'spec',       td = struct('height', 1.4,  'label', 'Freq (Hz)', 'order', 20);
    case 'hypnogram',  td = struct('height', 0.28, 'label', 'State',     'order', 10);
    case 'eventTicks', td = struct('height', 0.28, 'label', 'Events',    'order', 40);
    case 'stateStrip', td = struct('height', 0.5,  'label', 'State',     'order', 15);
    case 'raster',     td = struct('height', 1.2,  'label', 'Units',     'order', 60);
    case 'trace',      td = struct('height', 1.0,  'label', '',          'order', 50);
    case 'traces',     td = struct('height', 2.5,  'label', 'LFP',       'order', 50);
    otherwise,         td = struct('height', 1.0,  'label', '',          'order', 99);
end
end
