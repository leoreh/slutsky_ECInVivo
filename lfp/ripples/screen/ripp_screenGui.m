function hFig = ripp_screenGui(res, varargin)
% RIPP_SCREENGUI Compare methods' ripple events in one guiPath window.
%
%   hFig = RIPP_SCREENGUI(res, varargin)
%
%   SUMMARY:
%       Opens guiPath over one session with the ripple context (stacked-shank
%       LFP, filtered trace, EMG, unit raster, hypnogram, spectrogram) and every
%       method's detected events overlaid, each in its own colour: a tick strip
%       per method in the overview and colour-matched marks over the LFP. Uses
%       guiPath's multi-set support - each method is one event SET (a top tick
%       strip and a bottom overlay that share a name), so guiPath assigns the
%       colours and lists the methods under CURATE, and all methods sit side by
%       side in one window. Pick a method from CURATE to step its events
%       (arrows) with the others still overlaid; CURATE = None just steps the
%       window. Nothing is saved - the save target is neutralised, so
%       <basename>.ripp.mat is never touched.
%
%   INPUTS:
%       res      - <struct> ripp_screen output (needs .pk, .detect, .methods,
%                           .basepaths, .sbjID).
%       varargin - Parameter/Value:
%           'mouse'   - <num>  mouse index into res. {max disagreement}
%           'Visible' - <char> 'on' | 'off' for headless / scripted use. {'on'}
%
%   OUTPUT:
%       hFig - <handle> the guiPath figure.
%
%   DEPENDENCIES:
%       guiPath, guiPath_preset, guiPath_panel.
%
%   HISTORY:
%       260716 detection-review parameter screen.
%       260717 rewritten onto guiPath's multi-set events (one set per method).

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addParameter(p, 'mouse', [], @(x) isempty(x) || isscalar(x));
addParameter(p, 'Visible', 'on', @(x) any(strcmpi(char(x), {'on', 'off'})));
parse(p, varargin{:});
mouse = p.Results.mouse;
vis   = char(p.Results.Visible);

names = {res.methods.name};
nM    = numel(names);

%% ========================================================================
%  PICK THE MOUSE WHERE THE METHODS DISAGREE MOST
%  ========================================================================

if isempty(mouse)
    iCur = max(1, find(strcmp(names, 'default'), 1));
    iAlt = find(strcmp(names, 'fooof'), 1);
    if isempty(iAlt), iAlt = min(2, nM); end
    dis = zeros(1, numel(res.sbjID));
    for iMouse = 1:numel(res.sbjID)
        a = res.pk{iMouse, iCur};
        b = res.pk{iMouse, iAlt};
        dis(iMouse) = nUnmatched(a, b, 0.03) + nUnmatched(b, a, 0.03);
    end
    [~, mouse] = max(dis);
    fprintf('[SCREENGUI] most disagreement: mouse %d (%s), %d unmatched\n', ...
        mouse, res.sbjID{mouse}, dis(mouse));
end
basepath = res.basepaths{mouse};

%% ========================================================================
%  RIPPLE CONTEXT FROM THE PRESET, MINUS ITS SINGLE EVENTS PANEL
%  ========================================================================

[varMap, guiMap] = guiPath_preset('ripp', basepath);
if isfield(varMap, 'ripp'),       varMap = rmfield(varMap, 'ripp'); end
if isfield(guiMap.panels, 'evt'), guiMap.panels = rmfield(guiMap.panels, 'evt'); end

%% ========================================================================
%  ONE EVENT SET PER METHOD (own colour, top strip + bottom overlay)
%  ========================================================================
% Each method is a guiPath event set: a varMap entry (the events, inline) drawn
% by a top tick strip and a bottom overlay that share its var, so one coloured
% input shows in both regions. A method with no events is skipped.

for iMethod = 1:nM
    rp = res.detect{iMethod}(mouse).ripp;
    if isempty(rp) || ~any(rp.accepted)
        continue;
    end
    acc = rp.accepted;
    ev = struct('peakTime', rp.peakTime(acc), ...
        'times', rp.times(acc, :), ...
        'accepted', true(sum(acc), 1));
    fld = matlab.lang.makeValidName(names{iMethod});
    varMap.(fld) = var_recipe('value', 'data', ev);
    guiMap.panels.([fld, '_top']) = guiPath_panel('eventTicks', 'top', fld, ...
        'label', names{iMethod});
    guiMap.panels.([fld, '_bot']) = guiPath_panel('eventTicks', 'bottom', fld, ...
        'label', names{iMethod});
end

%% ========================================================================
%  LAUNCH; NO CURATION WRITES
%  ========================================================================

guiMap.name = 'Screen';
guiMap.save = '';               % neutralise the preset's 'ripp' save target
hFig = guiPath(basepath, 'varMap', varMap, 'guiMap', guiMap, 'Visible', vis);

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
for iVal = 1:numel(a)
    [~, iNear] = min(abs(b - a(iVal)));
    if abs(b(iNear) - a(iVal)) > tol
        n = n + 1;
    end
end
end
