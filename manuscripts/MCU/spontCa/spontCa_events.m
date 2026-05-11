function tbl = spontCa_events(tbl, fs, varargin)
% SPONTCA_EVENTS Runs single-trace detection on every row of a long-format
% SpontCa table and attaches per-event vectors plus per-cell scalars.
%
%   tbl = SPONTCA_EVENTS(tbl, fs, ...) runs SPONTCA_DETECT independently on
%   each row's trace. Cyto and mito are detected with separately tunable
%   parameters (via 'paramsCyto' and 'paramsMito'). Coupling between
%   compartments is NOT computed here - see SPONTCA_COUPLE.
%
%   COLUMNS ADDED:
%       start    (cell)    per-event start times (s)
%       stop     (cell)    per-event stop times (s)
%       amp      (cell)    per-event peak amp (dF/F)
%       dur      (cell)    per-event duration (s)
%       int      (cell)    per-event integral (dF/F * s)
%       noise    (n x 1)   robust noise std (from SPONTCA_DETECT)
%       nEvents  (n x 1)   event count
%       rate     (n x 1)   events / second
%       flux     (n x 1)   sum(amp) / recDur  (= rate * meanAmp)
%       fluxInt  (n x 1)   time-average of the trace
%       mapCyto  (cell)    [nCytoEv x nWin] this row's trace windowed
%                          around the cell's CYTO start times
%       mapMito  (cell)    [nMitoEv x nWin] this row's trace windowed
%                          around the cell's MITO start times
%
%   tbl.Properties.UserData carries 'tWin' (1 x nWin window time axis, s)
%   and 'mapWin' (the [pre post] window in seconds).
%
%   OPTIONAL (Name-Value):
%       'paramsCyto' - (cell) name-value args forwarded to SPONTCA_DETECT
%                      for cyto rows. Default {}.
%       'paramsMito' - (cell) name-value args forwarded to SPONTCA_DETECT
%                      for mito rows. Default {}.
%       'mapWin'     - (1x2) [pre post] window for ETA maps (s) {[-1, 5]}
%       'verbose'    - (log) print progress {true}.
%
%   See also: SPONTCA_LOAD, SPONTCA_DETECT, SPONTCA_COUPLE, SPONTCA_GUI

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'fs',  @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'paramsCyto', {}, @iscell);
addParameter(p, 'paramsMito', {}, @iscell);
addParameter(p, 'mapWin', [-1, 5], @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'verbose', true, @islogical);
parse(p, tbl, fs, varargin{:});
P = p.Results;

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
tbl.nEvents  = zeros(n, 1);
tbl.rate     = zeros(n, 1);
tbl.flux     = zeros(n, 1);
tbl.fluxInt  = zeros(n, 1);
tbl.mapCyto  = cell(n, 1);
tbl.mapMito  = cell(n, 1);

% Map window
winSamps = round(P.mapWin(1) * fs) : round(P.mapWin(2) * fs);
nWin     = length(winSamps);
tWin     = winSamps / fs;


%% ========================================================================
%  PER-ROW DETECTION
%  ========================================================================
% Detection is independent per row (per cell × compartment). Per-compartment
% params let the user tune cyto and mito separately.

for iRow = 1:n
    sig = tbl.trace(iRow, :);
    if all(isnan(sig))
        continue;
    end

    if tbl.compartment(iRow) == 'Cyto'
        detArgs = P.paramsCyto;
    else
        detArgs = P.paramsMito;
    end
    ev = spontCa_detect(sig, fs, detArgs{:});

    tbl.noise(iRow) = ev.noise;
    tbl.start{iRow} = ev.start;
    tbl.stop{iRow}  = ev.stop;
    tbl.amp{iRow}   = ev.amp;
    tbl.dur{iRow}   = ev.dur;
    tbl.int{iRow}   = ev.int;

    ampVec = ev.amp;
    valid  = ampVec > 0 & ~isnan(ampVec);
    nE     = sum(valid);
    tbl.nEvents(iRow) = nE;
    tbl.rate(iRow)    = nE / recDur;
    tbl.flux(iRow)    = sum(ampVec(valid)) / recDur;
    tbl.fluxInt(iRow) = mean(sig, 'omitnan');
end


%% ========================================================================
%  PER-CELL ETA MAPS (cyto-aligned AND mito-aligned)
%  ========================================================================

cells = unique(tbl.sbjID);

for iCell = 1:length(cells)
    sid = cells(iCell);
    iC  = find(tbl.sbjID == sid & tbl.compartment == 'Cyto');
    iM  = find(tbl.sbjID == sid & tbl.compartment == 'Mito');
    assert(isscalar(iC) && isscalar(iM), ...
        'spontCa_events: expect one row per compartment per cell');

    cyStarts = tbl.start{iC};
    miStarts = tbl.start{iM};

    % cyto-aligned maps for both rows
    if ~isempty(cyStarts)
        for iRow = [iC, iM]
            tbl.mapCyto{iRow} = ...
                buildMap(tbl.trace(iRow, :), cyStarts, dt, winSamps, nT);
        end
    end

    % mito-aligned maps for both rows
    if ~isempty(miStarts)
        for iRow = [iC, iM]
            tbl.mapMito{iRow} = ...
                buildMap(tbl.trace(iRow, :), miStarts, dt, winSamps, nT);
        end
    end
end

tbl.Properties.UserData.tWin   = tWin;
tbl.Properties.UserData.mapWin = P.mapWin;

if P.verbose
    fprintf('[spontCa_events] %d rows processed (independent detection)\n', n);
end

end     % EOF


%% ========================================================================
%  HELPER
%  ========================================================================

function rawMap = buildMap(sig, startTimes, dt, winSamps, nT)
% Vectorized event-triggered window extraction (RIPP_MAPS-style).
peakSamps = round(startTimes(:) / dt) + 1;
idxMat    = peakSamps + winSamps;
validMask = (idxMat >= 1) & (idxMat <= nT);
idxMat(~validMask) = 1;
rawMap = sig(idxMat);
rawMap(~validMask) = NaN;
end
