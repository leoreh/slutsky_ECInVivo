function [tblCell, tblEvent] = spontCa_finalize(tblCell, tblEvent, fs, varargin)
% SPONTCA_FINALIZE  Per-cell aggregates + ETA maps + mito-to-cyto coupling.
%
% Operates on the two-table representation:
%   tblCell  - one row per (cell, compartment); has trace + metadata.
%   tblEvent - one row per event; columns include sbjID, compartment,
%              start, stop, amp, dur, int.
%
% Output adds to tblCell: nEvents, rate, meanAmp, flux, fluxInt,
% mapCyto, mapMito, fracIndep.
%
% Output adds to tblEvent (mito rows only - cyto rows get NaN/false):
% cytoEvIdx, lag, cytoIndependent.
%
% Cyto stops in tblEvent are forced to start + 2*dt (dur = 2*dt, int = 0)
% because cyto decay times are not biologically real at fs=3 Hz.
%
% DESIGN RATIONALE - INDEPENDENT + COUPLING vs CYTO-TRIGGERED
% -----------------------------------------------------------
% The earlier "cyto-triggered" scheme generated one mito row per cyto
% event by definition, then declared "coupled" any window where the mito
% trace crossed threshold. This produced near-100% coupling rates and
% inflated mito event counts on flat traces. Independent detection lets
% each compartment speak for itself; the coupling step then asks: per
% mito event, was a cyto event "responsible"? The cytoIndependent
% fraction becomes the control.
%
% OPTIONAL (Name-Value):
%   'thrLag'  - lag above which a mito event is cyto-independent. {3 s}
%   'mapWin'  - [pre post] window for ETA maps (s). Default [-1, 5].
%   'verbose' - print summary. Default true.
%
% See also: SPONTCA_LOAD, SPONTCA_DETECT, SPONTCA_MANCUR, SPONTCA_GUI

%% ARGUMENTS

p = inputParser;
addRequired(p, 'tblCell',  @istable);
addRequired(p, 'tblEvent', @istable);
addRequired(p, 'fs',       @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'thrLag',  3,        @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'mapWin',  [-1, 5],  @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'verbose', true,     @islogical);
parse(p, tblCell, tblEvent, fs, varargin{:});
P = p.Results;

n      = height(tblCell);
nT     = size(tblCell.trace, 2);
recDur = nT / fs;
dt     = 1 / fs;


%% NORMALIZE CYTO STOPS
% Force cyto rows in tblEvent to placeholder stop/dur/int. Applies even
% if the user (or LLM) edited those values upstream - cyto decays aren't
% biologically real at fs=3.

isCyto = tblEvent.compartment == 'Cyto';
if any(isCyto)
    tblEvent.stop(isCyto) = tblEvent.start(isCyto) + 2 * dt;
    tblEvent.dur(isCyto)  = 2 * dt;
    tblEvent.int(isCyto)  = 0;
end


%% PER-ROW AGGREGATES

tblCell.nEvents = zeros(n, 1);
tblCell.rate    = zeros(n, 1);
tblCell.meanAmp = nan(n, 1);
tblCell.flux    = zeros(n, 1);
tblCell.fluxInt = zeros(n, 1);

for iRow = 1:n
    sig = tblCell.trace(iRow, :);
    tblCell.fluxInt(iRow) = mean(sig, 'omitnan');
    if all(isnan(sig)), continue; end

    mask = tblEvent.sbjID == tblCell.sbjID(iRow) & ...
           tblEvent.compartment == tblCell.compartment(iRow);
    amps = tblEvent.amp(mask);
    amps = amps(amps > 0 & ~isnan(amps));
    nE   = numel(amps);
    tblCell.nEvents(iRow) = nE;
    tblCell.rate(iRow)    = nE / recDur;
    if nE > 0
        tblCell.meanAmp(iRow) = mean(amps);
        tblCell.flux(iRow)    = sum(amps) / recDur;
    end
end


%% ETA MAPS

winSamps = round(P.mapWin(1) * fs) : round(P.mapWin(2) * fs);
tWin     = winSamps / fs;

tblCell.mapCyto = cell(n, 1);
tblCell.mapMito = cell(n, 1);

cells = unique(tblCell.sbjID, 'stable');
for iCell = 1:numel(cells)
    sid = cells(iCell);
    iC = find(tblCell.sbjID == sid & tblCell.compartment == 'Cyto');
    iM = find(tblCell.sbjID == sid & tblCell.compartment == 'Mito');
    if ~(isscalar(iC) && isscalar(iM))
        continue;
    end
    cyStarts = tblEvent.start(tblEvent.sbjID == sid & ...
        tblEvent.compartment == 'Cyto');
    miStarts = tblEvent.start(tblEvent.sbjID == sid & ...
        tblEvent.compartment == 'Mito');

    if ~isempty(cyStarts)
        for r = [iC, iM]
            tblCell.mapCyto{r} = ...
                buildMap(tblCell.trace(r, :), cyStarts, dt, winSamps, nT);
        end
    end
    if ~isempty(miStarts)
        for r = [iC, iM]
            tblCell.mapMito{r} = ...
                buildMap(tblCell.trace(r, :), miStarts, dt, winSamps, nT);
        end
    end
end

tblCell.Properties.UserData.tWin   = tWin;
tblCell.Properties.UserData.mapWin = P.mapWin;


%% COUPLING (mito -> preceding cyto)

nEv = height(tblEvent);
tblEvent.cytoEvIdx       = nan(nEv, 1);
tblEvent.lag             = inf(nEv, 1);
tblEvent.cytoIndependent = false(nEv, 1);
tblCell.fracIndep        = nan(n, 1);

for iCell = 1:numel(cells)
    sid = cells(iCell);
    cyMask = tblEvent.sbjID == sid & tblEvent.compartment == 'Cyto';
    miMask = tblEvent.sbjID == sid & tblEvent.compartment == 'Mito';
    if ~any(miMask), continue; end

    cyStarts = tblEvent.start(cyMask);
    miStarts = tblEvent.start(miMask);
    miRows   = find(miMask);

    cytoEvIdx = nan(numel(miStarts), 1);
    lag       = inf(numel(miStarts), 1);
    for m = 1:numel(miStarts)
        k = find(cyStarts <= miStarts(m), 1, 'last');
        if ~isempty(k)
            cytoEvIdx(m) = k;
            lag(m)       = miStarts(m) - cyStarts(k);
        end
    end
    indep = ~isfinite(lag) | lag > P.thrLag;

    tblEvent.cytoEvIdx(miRows)       = cytoEvIdx;
    tblEvent.lag(miRows)             = lag;
    tblEvent.cytoIndependent(miRows) = indep;

    iM = find(tblCell.sbjID == sid & tblCell.compartment == 'Mito');
    if isscalar(iM)
        tblCell.fracIndep(iM) = mean(indep);
    end
end


%% REPORT

if P.verbose
    iM = tblCell.compartment == 'Mito';
    fprintf(['[spontCa_finalize] %d cells | %d events total | ' ...
             'thrLag=%.1f s | median fracIndep = %.2f\n'], ...
        numel(cells), height(tblEvent), P.thrLag, ...
        median(tblCell.fracIndep(iM), 'omitnan'));
end

end     % EOF


function rawMap = buildMap(sig, startTimes, dt, winSamps, nT)
peakSamps = round(startTimes(:) / dt) + 1;
idxMat    = peakSamps + winSamps;
validMask = (idxMat >= 1) & (idxMat <= nT);
idxMat(~validMask) = 1;
rawMap = sig(idxMat);
rawMap(~validMask) = NaN;
end
