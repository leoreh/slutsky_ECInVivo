function [tblCell, tblEvent] = spontCa2_metrics(tblCell, tblEvent, fs, varargin)
% SPONTCA2_METRICS  Compute event + cell metrics from saved tables.
%
% [tblCell, tblEvent] = SPONTCA2_METRICS(tblCell, tblEvent, fs, ...)
% recomputes amp / dur / flux on every event from the per-cell traces,
% runs cross-compartment pairing, then aggregates to a cell-level
% summary. Designed to be re-callable on filtered subsets: the contract
% is tblEvent.sbjID is a subset of tblCell.sbjID; aggregates are
% computed only for cells present in tblCell.
%
% INPUTS
%   tblCell  - per-cell table with sbjID, compartment, genotype, trace.
%   tblEvent - long-format events: sbjID, compartment, start, stop
%              (and optionally genotype).
%   fs       - sampling rate (Hz).
%
% OPTIONAL (Name-Value)
%   'aggFcn'   - 'mean' (default) | 'median'. Per-cell aggregation of
%                event-level amp/dur/flux/pairFlux/pairAmp.
%   'winC2M'   - cyto -> mito response window (s). Default 2.
%   'maxLag'   - max cyto -> mito lag for trigger (s). Default 2.
%   'bslWin'   - rolling baseline window (s). Default 30. Mirrors the
%                detector for amp/flux recompute consistency.
%   'quantBsl' - rolling baseline percentile. Default 20.
%
% EVENT-LEVEL OUTPUT (overwrites or adds columns on tblEvent)
%   amp, dur, flux     - recomputed from trace using rolling 20th-pctile
%                        baseline. Cyto: dur=dt, flux=amp*dt, stop=start.
%                        Mito: dur=stop-start, flux=trapz(trace-bsl)*dt.
%   pairIdx, pairLag,
%   pairFlux, pairAmp  - cross-compartment pairing on the partner trace.
%   tf                 - pairFlux / flux (event-level transfer).
%   triggered          - categorical {'noPair','triggered'} for cyto rows.
%
% CELL-LEVEL OUTPUT (overwrites or adds columns on tblCell)
%   nEvents, rate      - count and event rate (Hz).
%   amp, dur, flux     - aggFcn over events.
%   pairFlux, pairAmp  - aggFcn over events.
%   fluxRate, ampRate  - sum(flux/amp) / recDur.
%   medAmp             - median(amp) regardless of aggFcn.
%   load               - sum(max(trace,0))*dt / recDur. Detection-free.
%
% PAIRING
%   cyto row: pairFlux/pairAmp = sum/max on mito trace over
%             [cyto.start, cyto.start + winC2M]. pairIdx points to the
%             first mito event whose start falls in the same window.
%   mito row: trigger = closest cyto event whose start is within
%             [mito.start - maxLag, mito.start]. pairFlux/pairAmp =
%             sum/max on cyto trace over [trigger.start, mito.stop].
%   Non-positive pairFlux/pairAmp/tf clamped to eps for log-friendly use.
%
% NOT computed here: per-cell transfer-function T = load_mito / load_cyto.
% That belongs to a dedicated TF analysis section in mcu_spontCa.
%
% See also: SPONTCA_DETECTWRAPPER, MCU_SPONTCA, SPONTCA_DETECT.


%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tblCell',  @istable);
addRequired(p, 'tblEvent', @istable);
addRequired(p, 'fs',       @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'aggFcn',   'mean', @(x) any(strcmpi(x, {'mean', 'median'})));
addParameter(p, 'winC2M',   2,      @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'maxLag',   2,      @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'bslWin',   30,     @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'quantBsl', 20,     @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 100);
parse(p, tblCell, tblEvent, fs, varargin{:});
P = p.Results;

switch lower(P.aggFcn)
    case 'mean',   aggH = @(x) mean(x, 'omitnan');
    case 'median', aggH = @(x) median(x, 'omitnan');
end

dt     = 1 / fs;
nSamps = size(tblCell.trace, 2);
recDur = nSamps / fs;
bslSmp = round(P.bslWin * fs);


%% ========================================================================
%  PREP: restrict events to cells present in tblCell
%  ========================================================================

isMatched = ismember(tblEvent.sbjID, tblCell.sbjID);
if any(~isMatched)
    fprintf('spontCa2_metrics: dropping %d events from cells not in tblCell\n', ...
        sum(~isMatched));
    tblEvent = tblEvent(isMatched, :);
end

nEv = height(tblEvent);


%% ========================================================================
%  EVENT-LEVEL: amp, dur, flux (recomputed from trace + baseline)
%  ========================================================================
% For each cell, build per-compartment rolling baselines once, then loop
% over its events. Cyto: stop forced to start; amp = trace - bsl at
% start; flux = amp*dt. Mito: amp = trace - bsl at start (which IS the
% peak by convention); flux = trapz(trace - bsl) over [start, stop] * dt.

tblEvent.amp  = nan(nEv, 1);
tblEvent.dur  = nan(nEv, 1);
tblEvent.flux = nan(nEv, 1);

cells = unique(tblEvent.sbjID);
for iCell = 1:numel(cells)
    sid = cells(iCell);
    rC = find(tblCell.sbjID == sid & tblCell.compartment == 'Cyto', 1);
    rM = find(tblCell.sbjID == sid & tblCell.compartment == 'Mito', 1);
    if isempty(rC) || isempty(rM), continue; end

    traceC = tblCell.trace(rC, :);
    traceM = tblCell.trace(rM, :);
    bslC   = rollingPercentile(traceC, bslSmp, P.quantBsl);
    bslM   = rollingPercentile(traceM, bslSmp, P.quantBsl);

    % --- Cyto events ---
    rowsC = find(tblEvent.sbjID == sid & tblEvent.compartment == 'Cyto');
    for k = 1:numel(rowsC)
        r  = rowsC(k);
        i0 = min(nSamps, max(1, round(tblEvent.start(r) * fs) + 1));
        tblEvent.stop(r) = tblEvent.start(r);   % enforce single-sample
        tblEvent.amp(r)  = traceC(i0) - bslC(i0);
        tblEvent.dur(r)  = dt;
        tblEvent.flux(r) = tblEvent.amp(r) * dt;
    end

    % --- Mito events ---
    rowsM = find(tblEvent.sbjID == sid & tblEvent.compartment == 'Mito');
    for k = 1:numel(rowsM)
        r  = rowsM(k);
        i0 = min(nSamps, max(1, round(tblEvent.start(r) * fs) + 1));
        i1 = min(nSamps, max(i0, round(tblEvent.stop(r) * fs) + 1));
        seg = traceM(i0:i1) - bslM(i0:i1);
        seg(isnan(seg)) = 0;
        tblEvent.amp(r)  = traceM(i0) - bslM(i0);
        tblEvent.dur(r)  = tblEvent.stop(r) - tblEvent.start(r);
        tblEvent.flux(r) = trapz(seg) * dt;
    end
end

% Defensive sanity print (never auto-drop here).
nanAmp = sum(isnan(tblEvent.amp));
if nanAmp > 0
    fprintf('spontCa2_metrics: %d events have NaN amp after recompute\n', nanAmp);
end


%% ========================================================================
%  EVENT-LEVEL: pairing + tf + triggered
%  ========================================================================

tblEvent.pairIdx  = nan(nEv, 1);
tblEvent.pairLag  = nan(nEv, 1);
tblEvent.pairFlux = nan(nEv, 1);
tblEvent.pairAmp  = nan(nEv, 1);
tblEvent.tf       = nan(nEv, 1);

for iCell = 1:numel(cells)
    sid = cells(iCell);
    rC = find(tblCell.sbjID == sid & tblCell.compartment == 'Cyto', 1);
    rM = find(tblCell.sbjID == sid & tblCell.compartment == 'Mito', 1);
    if isempty(rC) || isempty(rM), continue; end

    cytoTrace = tblCell.trace(rC, :);
    mitoTrace = tblCell.trace(rM, :);

    rowsC   = find(tblEvent.sbjID == sid & tblEvent.compartment == 'Cyto');
    rowsM   = find(tblEvent.sbjID == sid & tblEvent.compartment == 'Mito');
    startsC = tblEvent.start(rowsC);
    startsM = tblEvent.start(rowsM);
    stopsM  = tblEvent.stop(rowsM);

    % --- Cyto -> mito response window ---
    for k = 1:numel(rowsC)
        s  = startsC(k);
        i0 = max(1, round(s * fs) + 1);
        i1 = min(nSamps, round((s + P.winC2M) * fs) + 1);
        if i1 >= i0
            tblEvent.pairFlux(rowsC(k)) = sum(mitoTrace(i0:i1)) * dt;
            tblEvent.pairAmp(rowsC(k))  = max(mitoTrace(i0:i1));
        end
        j = find(startsM >= s & startsM <= s + P.winC2M, 1, 'first');
        if ~isempty(j)
            tblEvent.pairIdx(rowsC(k)) = rowsM(j);
            tblEvent.pairLag(rowsC(k)) = startsM(j) - s;
        end
    end

    % --- Mito -> trigger cyto, integrate cyto over [trigger.start, mito.stop] ---
    for k = 1:numel(rowsM)
        ms = startsM(k);
        j  = find(startsC <= ms & startsC >= ms - P.maxLag, 1, 'last');
        if isempty(j), continue; end
        trigStart = startsC(j);
        tblEvent.pairIdx(rowsM(k)) = rowsC(j);
        tblEvent.pairLag(rowsM(k)) = ms - trigStart;
        i0 = max(1, round(trigStart * fs) + 1);
        i1 = min(nSamps, round(stopsM(k) * fs) + 1);
        if i1 >= i0
            tblEvent.pairFlux(rowsM(k)) = sum(cytoTrace(i0:i1)) * dt;
            tblEvent.pairAmp(rowsM(k))  = max(cytoTrace(i0:i1));
        end
    end
end

% Per-event transfer ratio (units match: dF/F * s for both compartments).
tblEvent.tf = tblEvent.pairFlux ./ tblEvent.flux;

% Clamp non-positive entries to eps so log transforms survive.
tblEvent.pairFlux(tblEvent.pairFlux <= 0) = eps;
tblEvent.pairAmp(tblEvent.pairAmp   <= 0) = eps;
tblEvent.tf(tblEvent.tf <= 0)             = eps;

% Triggered categorical for cyto events.
isCytoEv = tblEvent.compartment == 'Cyto';
trig = NaN(nEv, 1);
trig(isCytoEv) = ~isnan(tblEvent.pairIdx(isCytoEv));
tblEvent.triggered = categorical(trig, [0 1], {'noPair', 'triggered'});


%% ========================================================================
%  CELL-LEVEL AGGREGATES
%  ========================================================================
% rate / fluxRate / ampRate are sum-over-time. amp / dur / flux /
% pairFlux / pairAmp use the configurable aggFcn. medAmp is always the
% median (kept for the compensation screen). load is detection-free.

nRows = height(tblCell);
tblCell.nEvents  = zeros(nRows, 1);
tblCell.rate     = zeros(nRows, 1);
tblCell.amp      = nan(nRows, 1);
tblCell.dur      = nan(nRows, 1);
tblCell.flux     = nan(nRows, 1);
tblCell.pairFlux = nan(nRows, 1);
tblCell.pairAmp  = nan(nRows, 1);
tblCell.fluxRate = zeros(nRows, 1);
tblCell.ampRate  = zeros(nRows, 1);
tblCell.medAmp   = nan(nRows, 1);
tblCell.load     = zeros(nRows, 1);

for iRow = 1:nRows
    mask = tblEvent.sbjID == tblCell.sbjID(iRow) & ...
           tblEvent.compartment == tblCell.compartment(iRow);
    n = sum(mask);
    tblCell.nEvents(iRow) = n;
    tblCell.rate(iRow)    = n / recDur;
    if n > 0
        tblCell.amp(iRow)      = aggH(tblEvent.amp(mask));
        tblCell.dur(iRow)      = aggH(tblEvent.dur(mask));
        tblCell.flux(iRow)     = aggH(tblEvent.flux(mask));
        tblCell.pairFlux(iRow) = aggH(tblEvent.pairFlux(mask));
        tblCell.pairAmp(iRow)  = aggH(tblEvent.pairAmp(mask));
        tblCell.fluxRate(iRow) = sum(tblEvent.flux(mask), 'omitnan') / recDur;
        tblCell.ampRate(iRow)  = sum(tblEvent.amp(mask),  'omitnan') / recDur;
        tblCell.medAmp(iRow)   = median(tblEvent.amp(mask), 'omitnan');
    end
    tr = tblCell.trace(iRow, :);
    tblCell.load(iRow) = sum(max(tr, 0)) * dt / recDur;
end

end     % SPONTCA2_METRICS


%% ========================================================================
%  HELPER: rolling percentile (centered window, NaN-tolerant, edge-shrunk)
%  ========================================================================
% Mirrors the helper inside spontCa_detect so amp/flux recomputation
% reproduces the detector's local baseline.

function out = rollingPercentile(x, winSmp, q)
nT = length(x);
out = nan(1, nT);
halfWin = floor(winSmp / 2);
for i = 1:nT
    lo = max(1, i - halfWin);
    hi = min(nT, i + halfWin);
    seg = x(lo:hi);
    seg = seg(~isnan(seg));
    if ~isempty(seg)
        out(i) = prctile(seg, q);
    end
end
end
