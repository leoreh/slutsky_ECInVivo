function pop = evt_pethPop(suMaps, muMaps, uType, tstamps)
% EVT_PETHPOP Population PETH (per event) from 3D spike maps, on demand.
%
%   pop = EVT_PETHPOP(suMaps, muMaps, uType, tstamps)
%
%   SUMMARY:
%       Reduces the 3D single-unit / MUA spike maps to per-event population
%       PETHs: one per cell type (RS, FS) plus pooled MUA, each the mean spike
%       train per event smoothed and z-scored against its own control
%       distribution (evt_pethNorm). Computed on demand from the maps stored in
%       the spks struct, so no population PETH is precomputed or saved.
%
%   INPUTS:
%       suMaps  - (Struct) Single-unit maps with .evt/.ctrl [N_units x N x nBins].
%       muMaps  - (Struct) Pooled-MUA maps with .evt/.ctrl [1 x N x nBins].
%       uType   - (Cat)    [N_units x 1] Unit types (RS/FS/...); [] to skip SU types.
%       tstamps - (Vec)    [1 x nBins] PETH time base (s).
%
%   OUTPUTS:
%       pop - (Struct) Per-event population PETH [N_events x nBins]:
%           .RS / .FS - Present only when that type has units.
%           .MU       - Pooled multi-unit (always).
%
%   DEPENDENCIES:
%       evt_pethNorm.
%
%   HISTORY:
%       Created: 05 Jul 2026 (replaces the population PETH block precomputed
%                into rippMaps.peth by ripp_wrapper).

pop = struct();

% Cell-type populations: mean spikes/unit per event, z-scored vs control
popTypes = {'RS', 'FS'};
for iType = 1:numel(popTypes)
    currType = popTypes{iType};
    if isempty(uType), continue; end
    idxType = uType == currType;
    if sum(idxType) == 0, continue; end

    nEvt  = size(suMaps.evt, 2);
    nCtrl = size(suMaps.ctrl, 2);
    nBins = size(suMaps.evt, 3);

    popEvt  = reshape(sum(suMaps.evt(idxType, :, :),  1), nEvt,  nBins) ./ sum(idxType);
    popCtrl = reshape(sum(suMaps.ctrl(idxType, :, :), 1), nCtrl, nBins) ./ sum(idxType);
    meanPopCtrl = mean(popCtrl, 1, 'omitnan');

    pop.(currType) = evt_pethNorm(popEvt, meanPopCtrl, tstamps);
end

% Pooled MUA
nEvt  = size(muMaps.evt, 2);
nCtrl = size(muMaps.ctrl, 2);
nBins = size(muMaps.evt, 3);

popEvt  = reshape(muMaps.evt,  nEvt,  nBins);
popCtrl = reshape(muMaps.ctrl, nCtrl, nBins);
meanPopCtrl = mean(popCtrl, 1, 'omitnan');
pop.MU = evt_pethNorm(popEvt, meanPopCtrl, tstamps);

end     % EOF
