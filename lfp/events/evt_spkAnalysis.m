function spks = evt_spkAnalysis(spkTimes, muTimes, evtTimes, ctrlTimes, peakTime, varargin)
% EVT_SPKANALYSIS Full spike-modulation analysis around events (one call).
%
%   spks = EVT_SPKANALYSIS(spkTimes, muTimes, evtTimes, ctrlTimes, peakTime, varargin)
%
%   SUMMARY:
%       Shared spike orchestration for both event pipelines (ripples, ED),
%       replacing the per-wrapper sequence of evt_spks + evt_spkPeth + kernel +
%       normalization. Combines, in one struct:
%       1. Per-unit scalar modulation stats + per-event population metrics
%          (.events), via evt_spks.
%       2. 3D peri-event spike maps for single units and pooled MUA, via
%          evt_spkPeth (.maps.su / .maps.mu, each .evt/.ctrl).
%       3. Per-unit normalized PETH (mean across events, smoothed + z-scored
%          against control), via evt_pethNorm (.peth).
%
%   INPUTS:
%       spkTimes  - (Cell) {N_units x 1} SU spike times (s).
%       muTimes   - (Cell) {1 x 1} Pooled MUA spike times (s).
%       evtTimes  - (Mat)  [N x 2] Event start/end times (s).
%       ctrlTimes - (Mat)  [N x 2] Matched control interval times (s).
%       peakTime  - (Vec)  [N x 1] Event peak times (s).
%       varargin  - Parameter/Value pairs:
%           'unitType' - (Cat) Unit types for per-type splits. (Default: []).
%           'mapDur'   - (Vec) PETH window [pre post] (s). (Default: [-0.05 0.05]).
%           'winFxd'   - (Num) Fixed half-window for scalar rates (s). (Default: 0.020).
%
%   OUTPUTS:
%       spks - (Struct) Consolidated spike results:
%           <per-unit scalar stats> - frEvt, frCtrl, cEvt, frZ, frMod, ... (evt_spks).
%           .events   - Per-event population metrics (frac/asym/com per type).
%           .peth     - [N_units x nBins] Per-unit normalized PETH.
%           .tstamps  - [1 x nBins] PETH time base (s).
%           .maps.su  - .evt/.ctrl 3D single-unit maps [N_units x N x nBins].
%           .maps.mu  - .evt/.ctrl 3D pooled-MUA maps [1 x N x nBins].
%
%   DEPENDENCIES:
%       evt_spks, evt_spkPeth, evt_pethNorm.
%
%   HISTORY:
%       Created: 05 Jul 2026 (absorbs the wrappers' spike orchestration).

% =========================================================================
%  ARGUMENTS
% =========================================================================
p = inputParser;
addRequired(p, 'spkTimes', @iscell);
addRequired(p, 'muTimes', @iscell);
addRequired(p, 'evtTimes', @isnumeric);
addRequired(p, 'ctrlTimes', @isnumeric);
addRequired(p, 'peakTime', @isnumeric);
addParameter(p, 'unitType', [], @(x) iscategorical(x) || iscell(x) || isempty(x));
addParameter(p, 'mapDur', [-0.05 0.05], @isnumeric);
addParameter(p, 'winFxd', 0.020, @(x) isnumeric(x) && isscalar(x) && x > 0);
parse(p, spkTimes, muTimes, evtTimes, ctrlTimes, peakTime, varargin{:});

unitType = p.Results.unitType;
mapDur   = p.Results.mapDur;
winFxd   = p.Results.winFxd;

% =========================================================================
%  PER-UNIT STATS + PER-EVENT POPULATION METRICS
% =========================================================================
spks = evt_spks(spkTimes, evtTimes, ctrlTimes, peakTime, ...
    'unitType', unitType, 'winFxd', winFxd, 'flgSave', false);

% =========================================================================
%  3D PETH MAPS (single units + pooled MUA)
% =========================================================================
suMaps = evt_spkPeth(spkTimes, peakTime, ctrlTimes, 'mapDur', mapDur, 'flgSave', false);
muMaps = evt_spkPeth(muTimes,  peakTime, ctrlTimes, 'mapDur', mapDur, 'flgSave', false);

spks.tstamps      = suMaps.tstamps;
spks.maps.su.evt  = suMaps.evt;
spks.maps.su.ctrl = suMaps.ctrl;
spks.maps.mu.evt  = muMaps.evt;
spks.maps.mu.ctrl = muMaps.ctrl;

% =========================================================================
%  PER-UNIT NORMALIZED PETH (mean across events, z-scored vs control)
% =========================================================================
nUnits = size(suMaps.evt, 1);
nBins  = size(suMaps.evt, 3);

% Mean across events -> [N_units x nBins] (reshape stays valid at N_units == 1)
meanPeth = reshape(mean(suMaps.evt,  2, 'omitnan'), nUnits, nBins);
ctrlPeth = reshape(mean(suMaps.ctrl, 2, 'omitnan'), nUnits, nBins);

spks.peth = evt_pethNorm(meanPeth, ctrlPeth, spks.tstamps);

end     % EOF
