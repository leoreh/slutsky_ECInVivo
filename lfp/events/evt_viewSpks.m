function evt_viewSpks(evtMaps, evtSpks, uType, state, varargin)
% EVT_VIEWSPKS Interactive tblGUI viewers for event maps + spike PETHs.
%
%   EVT_VIEWSPKS(evtMaps, evtSpks, uType, state, varargin)
%
%   SUMMARY:
%       Launches the exploratory tblGUI_xy viewers shared by both event
%       pipelines (ripples, ED), so the two reach parity from one code path:
%       1. LFP maps per event (from evtMaps), grouped by vigilance state.
%       2. Population PETH per event (RS/FS/MU), computed ON DEMAND from the 3D
%          spike maps (evt_pethPop) - never precomputed or saved.
%       3. Per-unit normalized PETH, grouped by unit type.
%       Every viewer is wrapped so a plotting hiccup cannot abort the pipeline.
%
%   INPUTS:
%       evtMaps  - (Struct) LFP maps with per-event fields + .tstamps.
%       evtSpks  - (Struct) Consolidated spike results (from evt_spkAnalysis),
%                           carrying .peth, .tstamps, .maps.su/.mu. [] to skip
%                           the spike viewers (e.g. no sorted spikes).
%       uType    - (Cat)    [N_units x 1] Unit types, or [].
%       state    - (Cat)    [N_events x 1] Per-event vigilance state, or [].
%       varargin - Parameter/Value pairs:
%           'mapYVar' - (Char) LFP map field to show first. (Default: 'z').
%
%   OUTPUTS:
%       None (launches figures).
%
%   DEPENDENCIES:
%       evt_pethPop, tblGUI_xy.
%
%   HISTORY:
%       Created: 05 Jul 2026 (replaces the inline viewer block in ripp_wrapper).

p = inputParser;
addParameter(p, 'mapYVar', 'z', @ischar);
parse(p, varargin{:});
mapYVar = p.Results.mapYVar;

% --- 1. LFP maps (per event) ---
try
    tblMap = table();
    mapFlds = setdiff(fieldnames(evtMaps), {'tstamps'});
    for iFld = 1:numel(mapFlds)
        tblMap.(mapFlds{iFld}) = evtMaps.(mapFlds{iFld});
    end
    if ~isempty(state) && height(tblMap) == numel(state)
        tblMap.states = state(:);
    end
    tblGUI_xy(evtMaps.tstamps, tblMap, 'yVar', mapYVar, 'grpVar', 'states');
catch ME
    warning('evt_viewSpks:lfpMaps', 'Skipped LFP-map viewer: %s', ME.message);
end

% Spike viewers require the consolidated spike maps
if isempty(evtSpks) || ~isfield(evtSpks, 'maps')
    return;
end

% --- 2. Population PETH (per event), on demand from the 3D ---
try
    pop = evt_pethPop(evtSpks.maps.su, evtSpks.maps.mu, uType, evtSpks.tstamps);
    tblPop = table();
    if isfield(pop, 'RS'), tblPop.pethRs = pop.RS; end
    if isfield(pop, 'FS'), tblPop.pethFs = pop.FS; end
    tblPop.pethMu = pop.MU;
    if ~isempty(state) && height(tblPop) == numel(state)
        tblPop.states = state(:);
    end
    tblGUI_xy(evtSpks.tstamps, tblPop, 'yVar', 'pethMu', 'grpVar', 'states');
catch ME
    warning('evt_viewSpks:popPeth', 'Skipped population-PETH viewer: %s', ME.message);
end

% --- 3. Per-unit normalized PETH, grouped by unit type ---
try
    tblUnit = table();
    if ~isempty(uType) && numel(uType) == size(evtSpks.peth, 1)
        tblUnit.unitType = uType(:);
    end
    tblUnit.peth = evtSpks.peth;
    tblGUI_xy(evtSpks.tstamps, tblUnit, 'yVar', 'peth', 'grpVar', 'unitType');
catch ME
    warning('evt_viewSpks:unitPeth', 'Skipped per-unit-PETH viewer: %s', ME.message);
end

end     % EOF
