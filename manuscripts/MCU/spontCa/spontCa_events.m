function tbl = spontCa_events(tbl, fs, varargin)
% SPONTCA_EVENTS Detects events on a long-format SpontCa table.
%
%   tbl = SPONTCA_EVENTS(tbl, fs, ...) adds per-event vectors (cell-array
%   columns) and per-cell scalars to the table returned by SPONTCA_LOAD.
%
%   MODES (name-value 'mode'):
%       'cytoTrigger' (default) - run SPONTCA_DETECT on the Cyto trace,
%           then build Mito events by windowing each cyto event:
%               window = [cyto_start_i, min(cyto_start_i + thrLag,
%                                            cyto_start_{i+1}, nT)]
%               mito_start = first sample inside the window where mito
%                            exceeds kThr * mito_noise. Uncoupled if no
%                            such sample exists.
%               mito_stop  = first sample after the in-window mito peak
%                            where the signal returns to thrBsl * peak,
%                            bounded by the right end of the window.
%           Guarantees one mito row entry per cyto event.
%       'independent' - run SPONTCA_DETECT on cyto and mito separately.
%           Used as a validation pipeline for the 'mito is cyto-triggered'
%           assumption (Atoms/MCU/MCU compensation model.md).
%
%   COLUMNS ADDED:
%       start    (cell)    per-event start times (s)
%       stop     (cell)    per-event stop times (s)
%       amp      (cell)    per-event peak amp (dF/F)
%       dur      (cell)    per-event duration (s)
%       int      (cell)    per-event integral (dF/F * s)
%       lag      (cell)    Mito rows: mito_start - cyto_start (s)
%       coupled  (cell)    Mito rows: logical, mito peak above noise
%       nEvents  (n x 1)   count of coupled events
%       rate     (n x 1)   events / second
%       flux     (n x 1)   sum(amp) / recDur  (matches rate * meanAmp)
%       fluxInt  (n x 1)   time-average of the trace
%       noise    (n x 1)   robust noise std
%       map      (cell)    [nCytoEv x nWin] event-triggered window of the
%                          row's own trace, aligned to each cyto start
%                          time (RIPP_MAPS-style). nWin shared across rows.
%
%   tbl.Properties.UserData carries 'tWin' (1 x nWin window time axis, s)
%   and 'mapWin' (the [pre post] window in seconds).
%
%   meanAmp / meanDur / meanInt are intentionally NOT stored - compute
%   on demand: cellfun(@mean, tbl.amp), etc.
%
%   OPTIONAL (Name-Value):
%       'mode'   - (char) 'cytoTrigger' (default) or 'independent'
%       'thrLag' - (num)  cyto-triggered search window (s) {3}
%       'thrBsl' - (num)  return-to fraction of peak for mito stop {0.1}
%       'mapWin' - (1x2)  [pre post] window for STA map (s) {[-1, 5]}
%       Pass-through to SPONTCA_DETECT: 'kThr', 'minAmp', 'minDur', 'minIEI'.
%       'verbose' - (log) print progress {true}.
%
%   See also: SPONTCA_LOAD, SPONTCA_DETECT, SPONTCA_GUI

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'fs',  @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'mode', 'cytoTrigger', ...
    @(x) any(strcmpi(x, {'cytoTrigger', 'independent'})));
addParameter(p, 'thrLag', 3,   @isnumeric);
addParameter(p, 'thrBsl', 0.1, @isnumeric);
addParameter(p, 'kThr',   3,    @isnumeric);
addParameter(p, 'minAmp', 0.02, @isnumeric);
addParameter(p, 'minDur', 0.4,  @isnumeric);
addParameter(p, 'minIEI', 0.4,  @isnumeric);
addParameter(p, 'mapWin', [-1, 5], @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'verbose', true, @islogical);
parse(p, tbl, fs, varargin{:});
P = p.Results;

detArgs = {'kThr', P.kThr, 'minAmp', P.minAmp, ...
    'minDur', P.minDur, 'minIEI', P.minIEI};

n      = height(tbl);
nT     = size(tbl.trace, 2);
recDur = nT / fs;
dt     = 1 / fs;


%% ========================================================================
%  ALLOCATE COLUMNS
%  ========================================================================

tbl.noise    = nan(n, 1);
tbl.start    = cell(n, 1);
tbl.stop     = cell(n, 1);
tbl.amp      = cell(n, 1);
tbl.dur      = cell(n, 1);
tbl.int      = cell(n, 1);
tbl.lag      = cell(n, 1);
tbl.coupled  = cell(n, 1);
tbl.nEvents  = zeros(n, 1);
tbl.rate     = zeros(n, 1);
tbl.flux     = zeros(n, 1);
tbl.fluxInt  = zeros(n, 1);
tbl.map      = cell(n, 1);

% Map window
winSamps = round(P.mapWin(1) * fs) : round(P.mapWin(2) * fs);
nWin     = length(winSamps);
tWin     = winSamps / fs;


%% ========================================================================
%  PER-CELL LOOP
%  ========================================================================

cells = unique(tbl.sbjID);

for iCell = 1:length(cells)

    sid = cells(iCell);
    iC  = find(tbl.sbjID == sid & tbl.compartment == 'Cyto');
    iM  = find(tbl.sbjID == sid & tbl.compartment == 'Mito');
    assert(isscalar(iC) && isscalar(iM), ...
        'spontCa_events: expect one row per compartment per cell');

    cyTrace = tbl.trace(iC, :);
    miTrace = tbl.trace(iM, :);
    if all(isnan(cyTrace)) || all(isnan(miTrace))
        continue;
    end

    % Cyto detection (always independent)
    evC = spontCa_detect(cyTrace, fs, detArgs{:});

    % Mito detection
    if strcmpi(P.mode, 'independent')
        evM = spontCa_detect(miTrace, fs, detArgs{:});
        lagVec  = nan(length(evM.start), 1);
        coupVec = true(length(evM.start), 1);
    else
        [evM, lagVec, coupVec] = ...
            cytoTriggeredMito(miTrace, evC.start, fs, P, dt, nT);
    end

    % Store cyto
    tbl.noise(iC)   = evC.noise;
    tbl.start{iC}   = evC.start;
    tbl.stop{iC}    = evC.stop;
    tbl.amp{iC}     = evC.amp;
    tbl.dur{iC}     = evC.dur;
    tbl.int{iC}     = evC.int;
    tbl.lag{iC}     = nan(length(evC.start), 1);
    tbl.coupled{iC} = false(length(evC.start), 1);

    % Store mito
    tbl.noise(iM)   = evM.noise;
    tbl.start{iM}   = evM.start;
    tbl.stop{iM}    = evM.stop;
    tbl.amp{iM}     = evM.amp;
    tbl.dur{iM}     = evM.dur;
    tbl.int{iM}     = evM.int;
    tbl.lag{iM}     = lagVec;
    tbl.coupled{iM} = coupVec;

    % Per-cell scalars for both rows
    for iRow = [iC, iM]
        ampVec  = tbl.amp{iRow};
        valid   = ampVec > 0 & ~isnan(ampVec);
        nE      = sum(valid);
        tbl.nEvents(iRow) = nE;
        tbl.rate(iRow)    = nE / recDur;
        tbl.flux(iRow)    = sum(ampVec(valid)) / recDur;
        tbl.fluxInt(iRow) = mean(tbl.trace(iRow, :), 'omitnan');
    end

    % --- STA map (RIPP_MAPS-style): windowed segments aligned to cyto
    %     start times, applied to BOTH cyto and mito traces of this cell.
    cytoStart = tbl.start{iC};
    if ~isempty(cytoStart)
        peakSamps = round(cytoStart / dt) + 1;
        idxMat = peakSamps + winSamps;
        validMask = (idxMat >= 1) & (idxMat <= nT);
        idxMat(~validMask) = 1;
        for iRow = [iC, iM]
            sig = tbl.trace(iRow, :);
            rawMap = sig(idxMat);
            rawMap(~validMask) = NaN;
            tbl.map{iRow} = rawMap;
        end
    end
end

% Stash window metadata for SPONTCA_GUI
tbl.Properties.UserData.tWin   = tWin;
tbl.Properties.UserData.mapWin = P.mapWin;

if P.verbose
    fprintf('[spontCa_events] mode=%s | %d cells processed\n', ...
        P.mode, length(cells));
end

end     % EOF


%% ========================================================================
%  HELPER (called twice: once per cytoTrigger pass + by independent mode
%  reusing the noise-floor logic from SPONTCA_DETECT directly)
%  ========================================================================

function [evM, lagVec, coupVec] = cytoTriggeredMito(miTrace, cStart, fs, P, dt, nT)
% Build mito events from cyto event windows. One row per cyto event.

% Borrow noise from a one-pass spontCa_detect call on the mito trace; the
% events list it returns is discarded - we only want the noise estimate.
miEv = spontCa_detect(miTrace, fs, ...
    'kThr', P.kThr, 'minAmp', P.minAmp, ...
    'minDur', P.minDur, 'minIEI', P.minIEI);
sNoise = miEv.noise;
thrMi  = max(P.kThr * sNoise, P.minAmp);

nC = length(cStart);
evM.start = nan(nC, 1);
evM.stop  = nan(nC, 1);
evM.amp   = nan(nC, 1);
evM.dur   = nan(nC, 1);
evM.int   = nan(nC, 1);
lagVec    = nan(nC, 1);
coupVec   = false(nC, 1);

cStartIdx = round(cStart / dt) + 1;
lagSmp    = round(P.thrLag / dt);

for iE = 1:nC
    iWinL = cStartIdx(iE);
    iWinR = min(iWinL + lagSmp, nT);
    if iE < nC
        iWinR = min(iWinR, cStartIdx(iE + 1) - 1);
    end
    if iWinR <= iWinL
        evM.start(iE) = (iWinL - 1) * dt;
        evM.stop(iE)  = evM.start(iE);
        evM.dur(iE)   = 0;
        evM.int(iE)   = 0;
        continue;
    end

    seg     = miTrace(iWinL:iWinR);
    segMask = seg > thrMi;
    if ~any(segMask)
        evM.start(iE) = (iWinL - 1) * dt;
        evM.stop(iE)  = evM.start(iE);
        evM.dur(iE)   = 0;
        evM.int(iE)   = 0;
        continue;
    end

    % first crossing of noise floor inside the window
    ascRel = find(segMask, 1, 'first');
    mStart = iWinL + ascRel - 1;

    % peak inside window
    segValid = seg;
    segValid(isnan(segValid)) = -inf;
    [pkVal, pkRel] = max(segValid);
    iPk = iWinL + pkRel - 1;

    % stop: first return to thrBsl * peak after peak, bounded by window
    thrReturn = P.thrBsl * pkVal;
    iR = iPk;
    while iR < iWinR && miTrace(iR + 1) > thrReturn
        iR = iR + 1;
    end
    mStop = iR;

    segE = miTrace(mStart:mStop);
    segE(isnan(segE)) = 0;

    evM.start(iE) = (mStart - 1) * dt;
    evM.stop(iE)  = (mStop  - 1) * dt;
    evM.amp(iE)   = pkVal;
    evM.dur(iE)   = evM.stop(iE) - evM.start(iE);
    evM.int(iE)   = trapz(segE) * dt;
    lagVec(iE)    = evM.start(iE) - cStart(iE);
    coupVec(iE)   = true;
end

evM.noise = sNoise;

end
