function hFig = ripp_screenGui(res, varargin)
% RIPP_SCREENGUI Inspect where detection methods disagree, in one guiPath window.
%
%   hFig = RIPP_SCREENGUI(res, varargin)
%
%   Opens guiPath over one session (default: the one where methods disagree
%   most) with the usual ripple context plus a raster panel that shows every
%   method's ripple peaks as its own row - all methods side by side in one
%   window. One method is the steppable target so the bottom window can walk
%   event to event. Saving is a no-op, so <basename>.ripp.mat is never touched.
%
%   The per-method peak times come straight from res.pk{mouse, method}, so
%   adding, dropping, or reordering a row is a one-line edit - no digging
%   through res.detect. guiPath itself is unchanged.
%
%   LIMITATION: the rows share one colour. guiPath's raster paints the whole
%   panel one colour and eventTicks only ever shows the single target, so
%   per-method colours would need a change inside the gui (drawRaster /
%   plot_raster), which this deliberately does not touch.
%
%   INPUTS:
%       res      - <struct> ripp_screen output (needs .pk, .methods, .sbjID).
%       varargin - Parameter/Value:
%           'mouse'  - <num>  mouse index into res. {max disagreement}
%           'target' - <char> method name used as the steppable target. {'narrow'}
%
%   OUTPUT:
%       hFig - <handle> the guiPath figure.
%
%   DEPENDENCIES:
%       guiPath, guiPath_presets.
%
%   HISTORY:
%       260716 detection-review parameter screen.

p = inputParser;
addParameter(p, 'mouse', [], @(x) isempty(x) || isscalar(x));
addParameter(p, 'target', 'narrow', @ischar);
parse(p, varargin{:});
mouse  = p.Results.mouse;
target = p.Results.target;

names = {res.methods.name};
nM    = numel(names);

%% ---- pick the mouse where 'current' and 'narrow' disagree most ----
if isempty(mouse)
    iCur = max(1, find(strcmp(names, 'current'), 1));
    iNar = find(strcmp(names, 'narrow'), 1);
    if isempty(iNar), iNar = min(2, nM); end
    dis = zeros(1, numel(res.sbjID));
    for iMouse = 1:numel(res.sbjID)
        a = res.pk{iMouse, iCur}; b = res.pk{iMouse, iNar};
        dis(iMouse) = nUnmatched(a, b, 0.03) + nUnmatched(b, a, 0.03);
    end
    [~, mouse] = max(dis);
    fprintf('[SCREENGUI] most disagreement: mouse %d (%s), %d unmatched\n', ...
        mouse, res.sbjID{mouse}, dis(mouse));
end
basepath = res.basepaths{mouse};

%% ---- context panels from the ripp preset, minus its curation target ----
cfgData = guiPath_presets('ripp');
if isfield(cfgData, 'evt'), cfgData = rmfield(cfgData, 'evt'); end

%% ---- one raster row per method (peaks straight from res.pk) ----
cfgData.methods = mkPanel('raster', 'bottom', '', 1.5, res.pk(mouse, :));
cfgData.methods.label = ['Methods top-down: ' strjoin(names, ' / ')];

%% ---- steppable target: one method's events (pre-loaded, not the .mat) ----
iTgt = find(strcmp(names, target), 1);
if isempty(iTgt), iTgt = min(2, nM); end
tgt = res.detect{iTgt}(mouse).ripp;
evtData = struct('peakTime', tgt.peakTime(:), 'times', tgt.times, ...
    'accepted', true(numel(tgt.peakTime), 1));
cfgData.evt = mkPanel('eventTicks', 'top', ['Stepper: ' names{iTgt}], 1, evtData);

%% ---- launch, save disabled ----
cfgGui = struct('name', 'Screen', 'file', '', 'mode', 'events', ...
    'win', 0.4, 'save', @(x) []);
hFig = guiPath(basepath, 'cfgData', cfgData, 'cfgGui', cfgGui);

end     % EOF


% =========================================================================
%  LOCALS
% =========================================================================
function n = nUnmatched(a, b, tol)
% count entries of a with no b within tol
n = 0;
if isempty(a), return; end
if isempty(b), n = numel(a); return; end
b = sort(b(:));
for i = 1:numel(a)
    [~, j] = min(abs(b - a(i)));
    if abs(b(j) - a(i)) > tol, n = n + 1; end
end
end

% -------------------------------------------------------------------------
function pnl = mkPanel(type, region, label, height, data)
% minimal panel struct guiPath_load accepts; .data preset -> not reloaded.
% {data} wraps the payload so struct() stays scalar even when data is a cell.
pnl = struct('type', type, 'region', region, 'src', [], 'name', label, ...
    'fs', NaN, 'height', height, 'clr', 'k', 'label', label, ...
    'order', 99, 'ylim', [], 'data', {data});
end
