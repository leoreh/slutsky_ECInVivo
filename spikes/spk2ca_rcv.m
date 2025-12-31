function [tbl, hFig] = spk2ca_rcv(v, varargin)
% SPK2CA_RCV Composite analysis of Spiking and Calcium recovery dynamics.
%
%   [tbl, hFig] = SPK2CA_RCV(v, ...) recalculates 1s binned metrics for
%   Rate, MinISI, CV, Cyto, and Mito. It then computes summary statistics
%   (Means, LogRatios) for Baseline and Steady-State windows and plots
%   the results.
%
%   INPUTS:
%       v           - (struct) Data structure containing .mea and .ca fields.
%       varargin    - (param/value) Optional parameters:
%                     'winBsl'  : (1x2) Baseline window (s) {[60, 4200]}
%                     'winSs'   : (1x2) Steady-State window (s) {[25200, 32395]}
%                     'winLim'  : (1x2) Overall analysis limits (s) {[0, 32400]}
%                     'flgPlot' : (log) Plot results {true}
%
%   OUTPUTS:
%       tbl         - (table) Unit table with traces (matrix cols) and
%                     stats (scalar cols).
%       hFig        - (handle) Figure handle of the plot.
%
%   See also: SPK2CA, TBLGUI_SCATHIST, MEA_TBL

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'v', @isstruct);
addParameter(p, 'winBsl', [1, 70] * 60, @isnumeric);
addParameter(p, 'winSs', [7 * 60, 9 * 60 - 5] * 60, @isnumeric);
addParameter(p, 'winLim', [0, 9 * 60] * 60, @isnumeric);
addParameter(p, 'flgPlot', true, @islogical);

parse(p, v, varargin{:});
winBsl = p.Results.winBsl;
winSs = p.Results.winSs;
winLim = p.Results.winLim;
flgPlot = p.Results.flgPlot;

%% ========================================================================
%  INITIALIZE
%  ========================================================================

% Flatten v if it's an array of files
if length(v) > 1
    v = v(:);
end
nFiles = length(v);

% Accumulate unit data
tblCells = cell(nFiles, 1);

fprintf('Processing %d Files...\n', nFiles);

%% ========================================================================
%  COMPUTE LOOP
%  ========================================================================

for iFile = 1

    % Extract Data
    if ~isfield(v(iFile), 'mea') || ~isfield(v(iFile).mea, 'spktimes')
        warning('File %d missing field "mea.spktimes". Skipping.', iFile);
        continue;
    end

    spktimes = v(iFile).mea.spktimes;
    nUnits = length(spktimes);

    % --- 1. Limit Spike Times ---
    spktimes = cellfun(@(x) x(x >= winLim(1) & x <= winLim(2)), ...
        spktimes, 'UniformOutput', false);

    % --- 2. Calculate Calcium (1s bins) ---
    ca = v(iFile).ca;

    % --- ALIGNMENT TO GLOBAL GRID ---
    % spk2ca may truncate the time vector based on data availability.
    % We must align to a global grid defined by winLim to ensure table stackability.

    % Define Global Grid (Bin Centers)
    edges = winLim(1) : 1 : winLim(2);
    tBins_ref = edges(1:end-1) + 0.5;
    nBins_ref = length(tBins_ref);

    % Initialize Matrices with NaNs (or Zeros)
    rateMat   = nan(nUnits, nBins_ref);
    isiMinMat = nan(nUnits, nBins_ref);
    cvMat     = nan(nUnits, nBins_ref);
    cytoMat   = nan(nUnits, nBins_ref);
    mitoMat   = nan(nUnits, nBins_ref);

    % Find mapping from ca.time to tBins_ref
    % We use nearest neighbor or direct matching since both are 1s bins
    [~, idxRef] = ismember(round(ca.time), round(tBins_ref));
    validMap = idxRef > 0;
    idxRef = idxRef(validMap);

    % Fill Ca Matrices
    if ~isempty(idxRef)
        cytoMat(:, idxRef) = ca.cyto(:, validMap);
        mitoMat(:, idxRef) = ca.mito(:, validMap);
    end

    % --- 3. Calculate Spiking Metrics (1s bins) ---
    for iUnit = 1:nUnits
        st = spktimes{iUnit};
        if isempty(st), continue; end

        % Discretize using GLOBAL edges
        binIdx = discretize(st, edges);
        valid = ~isnan(binIdx);
        binIdx = binIdx(valid);

        % Fill Rate (Count)
        % Reset valid bins to 0 first if we want to distinguish 0 from NaN (Missing)
        % Assuming that if we are processing this file, then for the duration overlapping
        % with winLim, we have data.
        % However, without precise start/stop times of the file, we can't be perfect.
        % We will just use accumarray which gives 0s for empty bins.

        % If we really want NaNs for "outside file", we keep init as NaN.
        % But we only update indices present in binIdx.
        % For rates, empty bins imply 0 Hz.
        % So let's initialize rateMat to 0 for the columns where we have Ca data?
        if ~isempty(idxRef)
            rateMat(iUnit, idxRef) = 0;
        end

        counts = accumarray(binIdx, 1, [nBins_ref 1]);

        % Assign counts to rateMat, respecting NaNs outside valid range?
        % If we strictly use idxRef to define "valid data range", we should only update those.
        % But spikes are filtered by winLim.
        % If a spike exists at time T, it implies T is valid.
        % Let's assign counts to all valid bins found.

        validBinIndices = find(counts > 0);
        rateMat(iUnit, validBinIndices) = counts(validBinIndices);

        % ISI Metrics
        for iBin = unique(binIdx)'
            s = st(binIdx == iBin);
            if length(s) > 1
                isi = diff(s);
                isiMinMat(iUnit, iBin) = min(isi);
                cvMat(iUnit, iBin) = std(isi) / mean(isi);
            end
        end
    end

    % Update timeVec reference for next step (Summary Stats)
    % Note: valid for ALL files now
    timeVec = tBins_ref;

    % --- 4. Package into Sub-Table ---
    t = table();
    t.FileIdx = repmat(iFile, nUnits, 1);
    t.UnitIdx = (1:nUnits)';

    % Helper extraction for Group/Name if present
    if isfield(v(iFile), 'fr') && isfield(v(iFile).fr, 'info')
        info = v(iFile).fr.info;
        if isfield(info, 'group')
            t.Group = repmat(string(info.group), nUnits, 1);
        end
    else
        % Fallback or check another location like user input structure
        t.Group = repmat("Unknown", nUnits, 1);
    end

    % Assign Traces
    t.Trace_Rate   = rateMat;
    t.Trace_MinISI = isiMinMat;
    t.Trace_CV     = cvMat;
    t.Trace_Cyto   = cytoMat;
    t.Trace_Mito   = mitoMat;

    tblCells{iFile} = t;
end

% Combine
tbl = vertcat(tblCells{:});

% Update timeVec to Reference
timeVec = tBins_ref;
nBins = nBins_ref;

% Return early logic moved to here
if isempty(tbl)
    warning('No units processed.');
    hFig = [];
    return;
end


%% ========================================================================
%  ANALYSIS: SUMMARY STATS
%  ========================================================================

% Define Indices for Windows
% Assumes all files share the same relative time mapping to the 1s bins
% 'timeVec' from the last iteration is used as reference.
% (Assuming winLim is consistent across files)

idxBsl = timeVec >= winBsl(1) & timeVec <= winBsl(2);
idxSs  = timeVec >= winSs(1)  & timeVec <= winSs(2);

traceVars = {'Trace_Rate', 'Trace_MinISI', 'Trace_CV', 'Trace_Cyto', 'Trace_Mito'};
statVars  = {'Rate', 'MinISI', 'CV', 'Cyto', 'Mito'};

for iVar = 1:length(traceVars)

    nameIn  = traceVars{iVar};
    nameOut = statVars{iVar};

    mat = tbl.(nameIn);

    % Means
    % Use 'omitnan' for ISI metrics which are often NaN in empty bins
    valBsl = mean(mat(:, idxBsl), 2, 'omitnan');
    valSs  = mean(mat(:, idxSs),  2, 'omitnan');

    % Assign Means
    tbl.([nameOut '_Bsl']) = valBsl;
    tbl.([nameOut '_Ss'])  = valSs;

    % Log Ratio
    % Add half minimum to avoid log(0).
    % For Rate/Cyto/Mito 0 is possible. For ISI metrics, they are > 0 if exist.

    % Find global non-zero min for epsilon
    allVal = [valBsl; valSs];
    posVal = allVal(allVal > 0);
    if isempty(posVal)
        epsVal = 1e-6;
    else
        epsVal = min(posVal) / 2;
    end

    ratio = log2((valSs + epsVal) ./ (valBsl + epsVal));

    tbl.([nameOut '_LogRatio']) = ratio;

end


%% ========================================================================
%  PLOTTING
%  ========================================================================

if flgPlot

    % Set Plotting Defaults
    % X: Log Ratio of Firing Rate
    % Y: Log Ratio of Mitochondrial Ca Accumulation

    hFig = tblGUI_scatHist(tbl, ...
        'xVar', 'Rate_LogRatio', ...
        'yVar', 'Mito_LogRatio', ...
        'grpVar', 'Group');

    set(hFig, 'Name', 'Spike-Ca Recovery Dynamics');

else
    hFig = [];
end


end     % EOF
