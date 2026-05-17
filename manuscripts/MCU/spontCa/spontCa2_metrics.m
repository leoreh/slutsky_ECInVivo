function [tblCell, tblEvent] = spontCa2_metrics(tblCell, tblEvent, fs, varargin)
% SPONTCA2_METRICS  Compute event + cell metrics from saved tables.
%
% [tblCell, tblEvent] = SPONTCA2_METRICS(tblCell, tblEvent, fs, ...)
% computes metrics in two passes:
%
%   Pass 1 (event-level scalars, idempotent under filtering): amp / dur /
%     flux / fluxOther / snr / traceSD / F0 / T. Derive from trace + per-
%     event start/stop alone. Run once on the post-detection table.
%
%   Pass 2 (event-set-dependent): pairing (pairIdx / pairLag / pairFlux /
%     pairAmp) and all cell-level aggregates. Re-run after every filter
%     step.
%
% INPUTS
%   tblCell  - per-cell table with sbjID, compartment, genotype, trace.
%   tblEvent - long-format events: sbjID, compartment, start, stop
%              (and optionally genotype).
%   fs       - sampling rate (Hz).
%
% OPTIONAL (Name-Value)
%   'mode'        - 'full' (default): runs Pass 1 + Pass 2.
%                   'cellOnly':       runs Pass 2 only. Auto-falls back to
%                                     'full' (with a printed note) if any
%                                     Pass 1 column is missing.
%   'aggFcn'      - 'mean' (default) | 'median'. Per-cell aggregation
%                   of event-level scalars (amp, dur, flux, pairFlux,
%                   pairAmp, fluxOther).
%   'winPair'     - 2-element vector [backward, forward] in seconds.
%                   Two roles: (a) fluxOther integration window for both
%                   compartments (see EVENT-LEVEL OUTPUT below); (b) the
%                   window that decides the paired logical column on
%                   tblEvent. Pairing itself is mutual nearest-neighbour
%                   with no window cut at pairing time. Default
%                   [0.34, 0.68] = 1 sample back and 2 samples forward
%                   at fs=3 Hz.
%
% EVENT-LEVEL OUTPUT (overwrites or adds columns on tblEvent)
%   amp, dur, flux    - read directly off trace. Cyto: dur=dt,
%                       flux=amp*dt, stop=start. Mito: dur=stop-start,
%                       flux=trapz(trace)*dt over [start, stop].
%   fluxOther         - detection-free integral of partner trace over a
%                       compartment-appropriate window:
%                         cyto row: [start - winPair(1), start + winPair(2)]
%                                   of mito trace.
%                         mito row: [start - winPair(2), start + winPair(1)]
%                                   of cyto trace (causally flipped - cyto
%                                   trigger is expected backward up to
%                                   winPair(2), with a small forward
%                                   window for jitter).
%                       Independent of partner-event detection.
%   traceSD           - per-cell SD of trace with event windows masked
%                       (start-1 through stop+1 samples). Computed against
%                       the full unfiltered event set in Pass 1 and locked.
%   snr               - amp / traceSD. Locked in Pass 1.
%   F0                - per-event lookup of tblCell.F0 when present
%                       (flgRaw=true path).
%   pairIdx, pairLag  - cross-compartment mutual nearest-neighbour
%                       partner. Sign convention is uniform across rows:
%                       pairLag = mito.start - cyto.start. Positive =
%                       mito follows cyto (causal). Negative = mito
%                       precedes cyto (jitter or atypical). pairIdx is
%                       NaN when no mutual partner exists (cell lacks
%                       partner-compartment events, or another event in
%                       the compartment was the mutual partner instead).
%   paired            - logical, true iff pairLag falls inside
%                       [-winPair(1), winPair(2)]. Decided at metrics
%                       time; spontCa_filter picks rows by this flag via
%                       flgPair. Mutual NN guarantees the paired cyto
%                       count equals the paired mito count.
%   pairFlux, pairAmp - partner event's event-level flux / amp (on the
%                       partner's own scale). NaN if no pair.
%   T                 - log transfer-function (output/input) per event:
%                         cyto row: log(fluxOther / flux)  = log(mito
%                                   response / cyto event)
%                         mito row: log(flux / fluxOther)  = log(mito
%                                   event / cyto context)
%                       Paired log ratio cancels event-level common
%                       variation (shared upstream drive).
%
% CELL-LEVEL OUTPUT (overwrites or adds columns on tblCell)
%   nEvents, rate                     - count and event rate (Hz).
%   amp, dur, flux                    - aggFcn over events.
%   pairFlux, pairAmp, fluxOther      - aggFcn over events.
%   fluxRate, ampRate                 - sum(flux/amp) / recDur.
%   load                              - sum(max(trace,0))*dt / recDur.
%                                       Detection-free.
%   T_cyto, T_mito                    - aggFcn over T for cyto events
%                                       and mito events respectively.
%                                       Attached to both rows of each
%                                       cell pair (identical values
%                                       across Cyto and Mito rows of the
%                                       same sbjID) for wide-format use.
%
% PAIRING
%   Mutual nearest-neighbour. Each event k in compartment A pairs with
%   event j in compartment B only when k is also the nearest A-event to
%   j. Non-mutual matches keep pairIdx NaN. Guarantees a strict 1:1
%   correspondence between paired cyto and paired mito events, even when
%   intra-compartment clustering brings two events within the winPair
%   span of the same partner. The paired flag (pairLag inside winPair)
%   is stored on tblEvent; spontCa_filter selects by it via flgPair.
%
%   Non-positive pairFlux/pairAmp/fluxOther clamped to eps for log-
%   friendly use.
%
% NOTE on baseline: this function does NOT subtract a rolling baseline
% from the trace. dF/F from spontCa_loadXls already includes the rolling
% 20th-percentile baseline; amp / flux / fluxOther are read directly off
% that trace.
%
% NOT computed here: per-cell transfer-function T = load_mito /
% load_cyto. That belongs to a dedicated TF analysis section in
% mcu_spontCa.
%
% See also: SPONTCA_DETECTWRAPPER, MCU_SPONTCA, SPONTCA_DETECT.


%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tblCell',  @istable);
addRequired(p, 'tblEvent', @istable);
addRequired(p, 'fs',       @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'mode',    'full', @(x) any(strcmpi(x, {'full', 'cellOnly'})));
addParameter(p, 'aggFcn',  'mean', @(x) any(strcmpi(x, {'mean', 'median'})));
addParameter(p, 'winPair', [0.34, 0.68], ...
    @(x) isnumeric(x) && numel(x) == 2 && all(x >= 0));
parse(p, tblCell, tblEvent, fs, varargin{:});
P = p.Results;

switch lower(P.aggFcn)
    case 'mean',   aggH = @(x) mean(x, 'omitnan');
    case 'median', aggH = @(x) median(x, 'omitnan');
end

dt     = 1 / fs;
nSamps = size(tblCell.trace, 2);
recDur = nSamps / fs;


%% ========================================================================
%  PREP: restrict events to cells present in tblCell + mode auto-fallback
%  ========================================================================

isMatched = ismember(tblEvent.sbjID, tblCell.sbjID);
if any(~isMatched)
    fprintf('spontCa2_metrics: dropping %d events from cells not in tblCell\n', ...
        sum(~isMatched));
    tblEvent = tblEvent(isMatched, :);
end

% Auto-fallback: cellOnly requires Pass 1 columns; otherwise run full.
needed = {'amp', 'dur', 'flux', 'fluxOther', 'snr'};
if strcmpi(P.mode, 'cellOnly') && ...
        ~all(ismember(needed, tblEvent.Properties.VariableNames))
    fprintf(['spontCa2_metrics: cellOnly requested but Pass 1 columns ', ...
        'missing on tblEvent; running full pass.\n']);
    P.mode = 'full';
end

nEv   = height(tblEvent);
cells = unique(tblEvent.sbjID);


%% ========================================================================
%  PASS 1: event-level scalars (amp, dur, flux, fluxOther, traceSD, snr, F0)
%  ========================================================================

if strcmpi(P.mode, 'full')

    tblEvent.amp       = nan(nEv, 1);
    tblEvent.dur       = nan(nEv, 1);
    tblEvent.flux      = nan(nEv, 1);
    tblEvent.fluxOther = nan(nEv, 1);
    tblEvent.snr       = nan(nEv, 1);

    nRows = height(tblCell);
    tblCell.traceSD = nan(nRows, 1);
    hasF0 = ismember('F0', tblCell.Properties.VariableNames);
    if hasF0
        tblEvent.F0 = nan(nEv, 1);
    end

    for iCell = 1:numel(cells)
        sid = cells(iCell);
        rC = find(tblCell.sbjID == sid & tblCell.compartment == 'Cyto', 1);
        rM = find(tblCell.sbjID == sid & tblCell.compartment == 'Mito', 1);
        if isempty(rC) || isempty(rM), continue; end

        traceC = tblCell.trace(rC, :);
        traceM = tblCell.trace(rM, :);

        rowsC = find(tblEvent.sbjID == sid & tblEvent.compartment == 'Cyto');
        rowsM = find(tblEvent.sbjID == sid & tblEvent.compartment == 'Mito');

        % --- Cyto event-level amp/dur/flux + fluxOther ---
        for k = 1:numel(rowsC)
            r  = rowsC(k);
            s  = tblEvent.start(r);
            i0 = min(nSamps, max(1, round(s * fs) + 1));
            tblEvent.stop(r) = s;
            tblEvent.amp(r)  = traceC(i0);
            tblEvent.dur(r)  = dt;
            tblEvent.flux(r) = tblEvent.amp(r) * dt;

            iF0 = max(1, round((s - P.winPair(1)) * fs) + 1);
            iF1 = min(nSamps, round((s + P.winPair(2)) * fs) + 1);
            if iF1 >= iF0
                tblEvent.fluxOther(r) = sum(traceM(iF0:iF1)) * dt;
            end
        end

        % --- Mito event-level amp/dur/flux + fluxOther ---
        for k = 1:numel(rowsM)
            r   = rowsM(k);
            ms  = tblEvent.start(r);
            me  = tblEvent.stop(r);
            i0  = min(nSamps, max(1, round(ms * fs) + 1));
            i1  = min(nSamps, max(i0, round(me * fs) + 1));
            seg = traceM(i0:i1);
            seg(isnan(seg)) = 0;
            tblEvent.amp(r)  = traceM(i0);
            tblEvent.dur(r)  = me - ms;
            tblEvent.flux(r) = trapz(seg) * dt;

            iF0 = max(1, round((ms - P.winPair(2)) * fs) + 1);
            iF1 = min(nSamps, round((ms + P.winPair(1)) * fs) + 1);
            if iF1 >= iF0
                tblEvent.fluxOther(r) = sum(traceC(iF0:iF1)) * dt;
            end
        end

        % --- Per-cell traceSD (event-masked) + per-event SNR + F0 ---
        for cmpIdx = 1:2
            if cmpIdx == 1
                tr   = traceC;
                rEv  = rowsC;
                rCll = rC;
            else
                tr   = traceM;
                rEv  = rowsM;
                rCll = rM;
            end
            sdMask = true(1, nSamps);
            for kE = 1:numel(rEv)
                iA = max(1, round(tblEvent.start(rEv(kE)) * fs) + 1 - 1);
                iB = min(nSamps, round(tblEvent.stop(rEv(kE))  * fs) + 1 + 1);
                sdMask(iA:iB) = false;
            end
            sdVal = std(tr(sdMask), 'omitnan');
            tblCell.traceSD(rCll) = sdVal;
            if sdVal > 0
                tblEvent.snr(rEv) = tblEvent.amp(rEv) / sdVal;
            end
            if hasF0
                tblEvent.F0(rEv) = tblCell.F0(rCll);
            end
        end
    end

    % Clamp non-positive fluxOther for log-friendly use.
    tblEvent.fluxOther(tblEvent.fluxOther <= 0) = eps;

    % Per-event log transfer-function T = log(output / input).
    %   Cyto row: mito response / cyto event   = log(fluxOther / flux)
    %   Mito row: mito event   / cyto context  = log(flux / fluxOther)
    tblEvent.T = nan(nEv, 1);
    isC = tblEvent.compartment == 'Cyto';
    isM = tblEvent.compartment == 'Mito';
    tblEvent.T(isC) = log(tblEvent.fluxOther(isC) ./ max(tblEvent.flux(isC), eps));
    tblEvent.T(isM) = log(max(tblEvent.flux(isM), eps) ./ tblEvent.fluxOther(isM));

    % Defensive sanity print.
    nanAmp = sum(isnan(tblEvent.amp));
    if nanAmp > 0
        fprintf('spontCa2_metrics: %d events have NaN amp after Pass 1\n', nanAmp);
    end
end


%% ========================================================================
%  PASS 2: pairing + cell aggregates
%  ========================================================================

tblEvent.pairIdx  = nan(nEv, 1);
tblEvent.pairLag  = nan(nEv, 1);
tblEvent.pairFlux = nan(nEv, 1);
tblEvent.pairAmp  = nan(nEv, 1);

for iCell = 1:numel(cells)
    sid = cells(iCell);
    rC = find(tblCell.sbjID == sid & tblCell.compartment == 'Cyto', 1);
    rM = find(tblCell.sbjID == sid & tblCell.compartment == 'Mito', 1);
    if isempty(rC) || isempty(rM), continue; end

    rowsC   = find(tblEvent.sbjID == sid & tblEvent.compartment == 'Cyto');
    rowsM   = find(tblEvent.sbjID == sid & tblEvent.compartment == 'Mito');
    startsC = tblEvent.start(rowsC);
    startsM = tblEvent.start(rowsM);
    if isempty(startsC) || isempty(startsM), continue; end

    % Mutual nearest-neighbour pairing: cyto k <-> mito j only when each
    % is the other's temporally closest event. Non-mutual matches leave
    % pairIdx NaN. Result: paired cyto count equals paired mito count
    % exactly, regardless of intra-compartment clustering. pairLag uses
    % the uniform sign convention mito.start - cyto.start.

    nearestM = zeros(numel(rowsC), 1);
    for k = 1:numel(rowsC)
        [~, nearestM(k)] = min(abs(startsM - startsC(k)));
    end
    nearestC = zeros(numel(rowsM), 1);
    for k = 1:numel(rowsM)
        [~, nearestC(k)] = min(abs(startsC - startsM(k)));
    end

    for k = 1:numel(rowsC)
        j = nearestM(k);
        if nearestC(j) == k
            mitoRow = rowsM(j);
            tblEvent.pairIdx(rowsC(k))  = mitoRow;
            tblEvent.pairLag(rowsC(k))  = startsM(j) - startsC(k);
            tblEvent.pairFlux(rowsC(k)) = tblEvent.flux(mitoRow);
            tblEvent.pairAmp(rowsC(k))  = tblEvent.amp(mitoRow);
        end
    end

    for k = 1:numel(rowsM)
        j = nearestC(k);
        if nearestM(j) == k
            cytoRow = rowsC(j);
            tblEvent.pairIdx(rowsM(k))  = cytoRow;
            tblEvent.pairLag(rowsM(k))  = startsM(k) - startsC(j);
            tblEvent.pairFlux(rowsM(k)) = tblEvent.flux(cytoRow);
            tblEvent.pairAmp(rowsM(k))  = tblEvent.amp(cytoRow);
        end
    end
end

% paired logical: pairLag inside winPair = [-winPair(1), +winPair(2)].
% Decided at metrics time using winPair; filters downstream pick rows
% by this flag instead of recomputing the window.
tblEvent.paired = tblEvent.pairLag >= -P.winPair(1) & ...
                  tblEvent.pairLag <= P.winPair(2);

tblEvent.pairFlux(tblEvent.pairFlux <= 0) = eps;
tblEvent.pairAmp(tblEvent.pairAmp   <= 0) = eps;


%% ========================================================================
%  CELL-LEVEL AGGREGATES
%  ========================================================================
% rate / fluxRate / ampRate are sum-over-time. amp / dur / flux /
% pairFlux / pairAmp / fluxOther use the configurable aggFcn. load is
% detection-free.

nRows = height(tblCell);
tblCell.nEvents   = zeros(nRows, 1);
tblCell.rate      = zeros(nRows, 1);
tblCell.amp       = nan(nRows, 1);
tblCell.dur       = nan(nRows, 1);
tblCell.flux      = nan(nRows, 1);
tblCell.pairFlux  = nan(nRows, 1);
tblCell.pairAmp   = nan(nRows, 1);
tblCell.fluxOther = nan(nRows, 1);
tblCell.fluxRate  = zeros(nRows, 1);
tblCell.ampRate   = zeros(nRows, 1);
tblCell.load      = zeros(nRows, 1);
tblCell.T_cyto    = nan(nRows, 1);
tblCell.T_mito    = nan(nRows, 1);

for iRow = 1:nRows
    sid  = tblCell.sbjID(iRow);
    mask = tblEvent.sbjID == sid & ...
           tblEvent.compartment == tblCell.compartment(iRow);
    n = sum(mask);
    tblCell.nEvents(iRow) = n;
    tblCell.rate(iRow)    = n / recDur;
    if n > 0
        tblCell.amp(iRow)       = aggH(tblEvent.amp(mask));
        tblCell.dur(iRow)       = aggH(tblEvent.dur(mask));
        tblCell.flux(iRow)      = aggH(tblEvent.flux(mask));
        tblCell.pairFlux(iRow)  = aggH(tblEvent.pairFlux(mask));
        tblCell.pairAmp(iRow)   = aggH(tblEvent.pairAmp(mask));
        tblCell.fluxOther(iRow) = aggH(tblEvent.fluxOther(mask));
        tblCell.fluxRate(iRow)  = sum(tblEvent.flux(mask), 'omitnan') / recDur;
        tblCell.ampRate(iRow)   = sum(tblEvent.amp(mask),  'omitnan') / recDur;
    end

    % T_cyto / T_mito: cell-mean T over each compartment's events,
    % attached to BOTH rows of the cell (identical across Cyto and Mito
    % rows of the same sbjID) so a single tblCell row holds the pair.
    maskC = tblEvent.sbjID == sid & tblEvent.compartment == 'Cyto';
    maskM = tblEvent.sbjID == sid & tblEvent.compartment == 'Mito';
    if any(maskC), tblCell.T_cyto(iRow) = aggH(tblEvent.T(maskC)); end
    if any(maskM), tblCell.T_mito(iRow) = aggH(tblEvent.T(maskM)); end

    tr = tblCell.trace(iRow, :);
    tblCell.load(iRow) = sum(max(tr, 0)) * dt / recDur;
end

end     % SPONTCA2_METRICS
