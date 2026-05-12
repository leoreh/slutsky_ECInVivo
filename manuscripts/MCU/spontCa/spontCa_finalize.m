function tbl = spontCa_finalize(tbl, fs, varargin)
% SPONTCA_FINALIZE Overlays curated events, computes per-row aggregates,
% per-cell ETA maps, and post-hoc mito-to-cyto coupling.
%
% Replaces the earlier SPONTCA_EVENTS + SPONTCA_COUPLE pair. The detection
% step (spontCa_detect per row) now lives inline in mcu_spontCa.m; this
% function consumes the resulting per-row event columns and produces
% everything needed for figures and QC.
%
% PIPELINE STEPS (in order):
%   1. CURATION OVERLAY. For each row, looks for
%      <curatedDir>/<sbjID>_<compartment>.mat (written by spontCa_manCur).
%      If present, overwrites tbl.start / .stop / .amp / .dur / .int for
%      that row with the curated values.
%   2. NORMALIZE CYTO STOPS. Cyto decay times are not biologically real
%      at fs=3 (subsequent events contaminate them). Cyto stop is forced
%      to start + 2 samples and dur/int are zeroed. Mito unchanged.
%   3. PER-ROW AGGREGATES. nEvents, rate, flux (sum(amp)/recDur), fluxInt
%      (time-mean of the trace).
%   3. ETA MAPS. Per cell, mapCyto (this row's trace windowed on the
%      cell's CYTO starts) and mapMito (windowed on the cell's MITO
%      starts). Stored on BOTH the cyto row and the mito row of the cell
%      so the QC viewer can overlay them.
%   4. COUPLING. Per mito event, finds the closest preceding cyto event.
%      Stores cytoEvIdx, lag, cytoIndependent (boolean), and per-cell
%      fracIndep (fraction of cyto-independent mito events).
%
% DESIGN RATIONALE - INDEPENDENT + COUPLING vs CYTO-TRIGGERED
% -----------------------------------------------------------
% The earlier "cyto-triggered" scheme generated one mito row per cyto
% event by definition, then declared "coupled" any window where the mito
% trace crossed threshold. This produced near-100% coupling rates and
% inflated mito event counts on flat traces (the window itself created
% the event). It also could not measure how often mito fires WITHOUT a
% cyto trigger - exactly the quantity needed to validate the model that
% mito is predominantly cyto-driven (Atoms/MCU/MCU compensation model.md).
%
% Independent detection lets each compartment speak for itself. The
% coupling step then asks the real question: per mito event, was a cyto
% event "responsible"? The cytoIndependent fraction becomes the control.
%
% OPEN CONCERN - COMPOUND EVENTS
% ------------------------------
% A single mito event may continue across several cyto triggers. Under
% cyto-triggered detection this loses per-trigger attribution. Two
% complementary safeguards:
%   (1) Event side: a fresh cyto trigger on a still-decaying mito should
%       produce a visible bump. Good mito detection (currently amp-only;
%       kinetics-based as future work) segments that bump as a new event.
%   (2) Population side: the cyto-aligned ETA averages mito amplitude
%       around every cyto onset, so even sub-peaks that detection misses
%       still contribute. ETA is therefore detection-quality-independent.
% These together make the cytoIndependent quantification credible without
% needing event detection to be perfect.
%
% COLUMNS ADDED:
%   nEvents          (n x 1)   event count per row
%   rate             (n x 1)   events / second
%   flux             (n x 1)   sum(amp) / recDur
%   fluxInt          (n x 1)   time-average of the trace
%   mapCyto          (cell)    [nCytoEv x nWin] this row's trace windowed
%                              on the cell's CYTO starts
%   mapMito          (cell)    [nMitoEv x nWin] this row's trace windowed
%                              on the cell's MITO starts
%   cytoEvIdx        (cell)    mito rows only: 1-based cyto event index
%                              (NaN if no preceding cyto event)
%   lag              (cell)    mito rows only: t_m - t_c (s); Inf if no
%                              preceding cyto event
%   cytoIndependent  (cell)    mito rows only: logical
%   fracIndep        (n x 1)   mito rows: fraction of cyto-independent
%                              events; cyto rows: NaN
%
% tbl.Properties.UserData carries 'tWin' (1 x nWin window time axis, s)
% and 'mapWin' (the [pre post] window in seconds).
%
% OPTIONAL (Name-Value):
%   'curatedDir' - (char) directory holding <sbjID>_<compartment>.mat
%                  files. Default: spontCa_curated/ alongside this file.
%   'thrLag'     - (num)  lag above which a mito event is cyto-independent
%                  {3}.
%   'mapWin'     - (1x2)  [pre post] window for ETA maps (s) {[-1, 5]}.
%   'verbose'    - (log)  print progress {true}.
%
% See also: SPONTCA_LOAD, SPONTCA_DETECT, SPONTCA_MANCUR, SPONTCA_GUI

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'fs',  @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'curatedDir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'thrLag', 3,    @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'mapWin', [-1, 5], @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'verbose', true, @islogical);
parse(p, tbl, fs, varargin{:});
P = p.Results;

if isempty(P.curatedDir)
    P.curatedDir = fullfile(fileparts(mfilename('fullpath')), ...
        'spontCa_curated');
end
P.curatedDir = char(P.curatedDir);

n      = height(tbl);
nT     = size(tbl.trace, 2);
recDur = nT / fs;
dt     = 1 / fs;


%% ========================================================================
%  CURATION OVERLAY
%  ========================================================================

nCurated = 0;
if isfolder(P.curatedDir)
    for iRow = 1:n
        sid = char(tbl.sbjID(iRow));
        cmp = char(tbl.compartment(iRow));
        fCur = fullfile(P.curatedDir, sprintf('%s_%s.mat', sid, cmp));
        if isfile(fCur)
            S = load(fCur);
            if isfield(S, 'cur')
                cur = S.cur;
                tbl.start{iRow} = cur.start(:);
                tbl.stop{iRow}  = cur.stop(:);
                tbl.amp{iRow}   = cur.amp(:);
                tbl.dur{iRow}   = cur.dur(:);
                tbl.int{iRow}   = cur.int(:);
                nCurated = nCurated + 1;
            end
        end
    end
end


%% ========================================================================
%  NORMALIZE CYTO STOPS
%  ========================================================================
% At fs=3 Hz with overlapping events, cyto decay times are contaminated
% by subsequent activity. The captured stop is not biologically real, so
% per-event dur/int aren't either. Force cyto stop = start + 2 samples
% (so the cell-array shapes stay consistent and downstream code that
% indexes them doesn't crash); zero out dur/int as a "do not interpret"
% marker. Applied AFTER curation overlay so old curation files with
% non-placeholder cyto stops get normalised the same way. Mito unchanged.

for iRow = 1:n
    if tbl.compartment(iRow) == 'Cyto'
        s = tbl.start{iRow};
        if ~isempty(s)
            tbl.stop{iRow} = s + 2 * dt;
            tbl.dur{iRow}  = repmat(2 * dt, numel(s), 1);
            tbl.int{iRow}  = zeros(numel(s), 1);
        end
    end
end


%% ========================================================================
%  PER-ROW AGGREGATES
%  ========================================================================

tbl.nEvents = zeros(n, 1);
tbl.rate    = zeros(n, 1);
tbl.flux    = zeros(n, 1);
tbl.fluxInt = zeros(n, 1);

for iRow = 1:n
    sig = tbl.trace(iRow, :);
    if all(isnan(sig))
        continue;
    end
    ampVec = tbl.amp{iRow};
    valid  = ampVec > 0 & ~isnan(ampVec);
    nE     = sum(valid);
    tbl.nEvents(iRow) = nE;
    tbl.rate(iRow)    = nE / recDur;
    tbl.flux(iRow)    = sum(ampVec(valid)) / recDur;
    tbl.fluxInt(iRow) = mean(sig, 'omitnan');
end


%% ========================================================================
%  ETA MAPS
%  ========================================================================

winSamps = round(P.mapWin(1) * fs) : round(P.mapWin(2) * fs);
tWin     = winSamps / fs;

tbl.mapCyto = cell(n, 1);
tbl.mapMito = cell(n, 1);

cells = unique(tbl.sbjID);
for iCell = 1:length(cells)
    sid = cells(iCell);
    iC  = find(tbl.sbjID == sid & tbl.compartment == 'Cyto');
    iM  = find(tbl.sbjID == sid & tbl.compartment == 'Mito');
    assert(isscalar(iC) && isscalar(iM), ...
        'spontCa_finalize: expect one row per compartment per cell');

    cyStarts = tbl.start{iC};
    miStarts = tbl.start{iM};

    if ~isempty(cyStarts)
        for iRow = [iC, iM]
            tbl.mapCyto{iRow} = ...
                buildMap(tbl.trace(iRow, :), cyStarts, dt, winSamps, nT);
        end
    end
    if ~isempty(miStarts)
        for iRow = [iC, iM]
            tbl.mapMito{iRow} = ...
                buildMap(tbl.trace(iRow, :), miStarts, dt, winSamps, nT);
        end
    end
end

tbl.Properties.UserData.tWin   = tWin;
tbl.Properties.UserData.mapWin = P.mapWin;


%% ========================================================================
%  COUPLING
%  ========================================================================

tbl.cytoEvIdx       = cell(n, 1);
tbl.lag             = cell(n, 1);
tbl.cytoIndependent = cell(n, 1);
tbl.fracIndep       = nan(n, 1);

for iCell = 1:length(cells)
    sid = cells(iCell);
    iC  = find(tbl.sbjID == sid & tbl.compartment == 'Cyto');
    iM  = find(tbl.sbjID == sid & tbl.compartment == 'Mito');

    cyStarts = tbl.start{iC};
    miStarts = tbl.start{iM};

    nM = length(miStarts);
    if nM == 0
        tbl.cytoEvIdx{iM}       = zeros(0, 1);
        tbl.lag{iM}             = zeros(0, 1);
        tbl.cytoIndependent{iM} = false(0, 1);
        tbl.fracIndep(iM)       = NaN;
        continue;
    end

    cytoEvIdx = nan(nM, 1);
    lag       = inf(nM, 1);
    for m = 1:nM
        k = find(cyStarts <= miStarts(m), 1, 'last');
        if ~isempty(k)
            cytoEvIdx(m) = k;
            lag(m)       = miStarts(m) - cyStarts(k);
        end
    end
    cytoIndependent = ~isfinite(lag) | lag > P.thrLag;

    tbl.cytoEvIdx{iM}       = cytoEvIdx;
    tbl.lag{iM}             = lag;
    tbl.cytoIndependent{iM} = cytoIndependent;
    tbl.fracIndep(iM)       = mean(cytoIndependent);
end


%% ========================================================================
%  REPORT
%  ========================================================================

if P.verbose
    iM = tbl.compartment == 'Mito';
    fprintf(['[spontCa_finalize] %d/%d rows used curated events | ' ...
             'thrLag=%.1f s | median fracIndep = %.2f\n'], ...
        nCurated, n, P.thrLag, median(tbl.fracIndep(iM), 'omitnan'));
end

end     % EOF


%% ========================================================================
%  HELPER
%  ========================================================================

function rawMap = buildMap(sig, startTimes, dt, winSamps, nT)
% Vectorized event-triggered window extraction.
peakSamps = round(startTimes(:) / dt) + 1;
idxMat    = peakSamps + winSamps;
validMask = (idxMat >= 1) & (idxMat <= nT);
idxMat(~validMask) = 1;
rawMap = sig(idxMat);
rawMap(~validMask) = NaN;
end
