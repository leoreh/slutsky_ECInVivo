function tbl = spontCa_couple(tbl, varargin)
% SPONTCA_COUPLE Post-hoc coupling of mito events to preceding cyto events.
%
% For each mito event with start t_m, finds the closest cyto event with
% start t_c <= t_m. Stores cytoEvIdx (the cyto event index, 1-based) and
% lag = t_m - t_c. If lag > thrLag, OR no cyto event precedes the mito
% event in the recording, the mito event is flagged cytoIndependent.
% Per-cell scalar fracIndep stores the fraction of cyto-independent events.
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
% COLUMNS ADDED (only meaningful on Mito rows; cyto rows hold empty/NaN):
%       cytoEvIdx        (cell)   nMitoEv x 1, 1-based cyto event index
%                                 (NaN if no preceding cyto event)
%       lag              (cell)   nMitoEv x 1, t_m - t_c (s); Inf if no
%                                 preceding cyto event
%       cytoIndependent  (cell)   nMitoEv x 1 logical
%       fracIndep        (n x 1)  per-cell fraction of cytoIndependent
%                                 (mito rows only; NaN on cyto rows)
%
% PARAMETERS (Name-Value):
%       'thrLag'  (3 s)  lag above which a mito event is cyto-independent
%       'verbose' (true) print summary
%
%   See also: SPONTCA_LOAD, SPONTCA_DETECT, SPONTCA_EVENTS, SPONTCA_GUI

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addParameter(p, 'thrLag', 3, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'verbose', true, @islogical);
parse(p, tbl, varargin{:});
P = p.Results;

n = height(tbl);


%% ========================================================================
%  ALLOCATE COLUMNS
%  ========================================================================

tbl.cytoEvIdx       = cell(n, 1);
tbl.lag             = cell(n, 1);
tbl.cytoIndependent = cell(n, 1);
tbl.fracIndep       = nan(n, 1);


%% ========================================================================
%  PER-CELL COUPLING
%  ========================================================================

cells = unique(tbl.sbjID);

for iCell = 1:length(cells)
    sid = cells(iCell);
    iC  = find(tbl.sbjID == sid & tbl.compartment == 'Cyto');
    iM  = find(tbl.sbjID == sid & tbl.compartment == 'Mito');
    assert(isscalar(iC) && isscalar(iM), ...
        'spontCa_couple: expect one row per compartment per cell');

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

if P.verbose
    iM = tbl.compartment == 'Mito';
    fprintf('[spontCa_couple] thrLag=%.1f s | median fracIndep = %.2f\n', ...
        P.thrLag, median(tbl.fracIndep(iM), 'omitnan'));
end

end     % EOF
