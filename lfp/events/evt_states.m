function [stateIdx, evtStates] = evt_states(evtTimes, peakTimes, boutTimes, varargin)
% EVT_STATES Classifies events by vigilance state and computes bout statistics.
%
%   [stateIdx, evtStates] = EVT_STATES(evtTimes, peakTimes, boutTimes, varargin)
%
%   SUMMARY:
%       1. Assigns a vigilance state to each event based on its peak time.
%       2. Calculates Event Rate (Hz) and Density (duration/duration) for
%          every individual sleep bout.
%       3. Aggregates this into a summary table.
%
%   INPUTS:
%       evtTimes   - (Mat)  [N x 2] Event start/end times [s].
%       peakTimes   - (Vec)  [N x 1] Event peak times [s].
%       boutTimes   - (Cell) {N_states x 1} Cell array of [Start End] matrices.
%                            Typically {Wake, NREM, REM, ...}.
%       varargin    - Parameter/Value pairs:
%           'basepath' - (Char) Save location. (Default: pwd).
%           'flgPlot'  - (Log)  Generate summary figures? (Default: true).
%           'flgSave'  - (Log)  Save .evtStates.mat? (Default: true).
%           'accepted' - (Log)  [N x 1] mask; Rate/Density count only accepted
%                              events, while stateIdx still labels all. {all}
%           'basename' - (Char) File stem for the saved table. {folder name}
%
%   OUTPUTS:
%       stateIdx    - (Cat)   [N_events x 1] Categorical array of states.
%       evtStates  - (Table) Summary table with rows per bout:
%                             [Rate, Density, Duration, State, Start, End].
%
%   DEPENDENCIES:
%       intervals, as_loadConfig.
%
%   HISTORY:
%       Updated: 23 Jan 2026

% =========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'evtTimes', @isnumeric);
addRequired(p, 'peakTimes', @isnumeric);
addRequired(p, 'boutTimes', @(x) iscell(x) || isempty(x));
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgSave', true, @islogical);
addParameter(p, 'accepted', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
addParameter(p, 'name', 'evt', @ischar);
addParameter(p, 'lbl', 'Event', @ischar);
addParameter(p, 'basename', '', @ischar);
parse(p, evtTimes, peakTimes, boutTimes, varargin{:});

evtTimes = p.Results.evtTimes;
peakTimes = p.Results.peakTimes;
boutTimes = p.Results.boutTimes;
if isempty(boutTimes), boutTimes = {}; end   % [] -> {} : no states, all events undefined
basepath = p.Results.basepath;
flgPlot = p.Results.flgPlot;
flgSave = p.Results.flgSave;
name = p.Results.name;
lbl = p.Results.lbl;

% =========================================================================
%  PREP
%  ========================================================================
% basename normally follows the folder, but several recordings can share one
% directory (the EA cohort), in which case the caller passes it explicitly -
% otherwise every session in the folder would overwrite one <folder>.States file
basename = p.Results.basename;
if isempty(basename)
    [~, basename] = fileparts(basepath);
end
savefile = fullfile(basepath, [basename, '.', name, 'States.mat']);

nStates = length(boutTimes);
nEvents = length(peakTimes);

% Acceptance mask: Rate/Density count only accepted events (default all)
accepted = p.Results.accepted;
if isempty(accepted)
    accepted = true(nEvents, 1);
else
    accepted = logical(accepted(:));
end

% Initialize Categorical Index
stateIdx = categorical(nan(nEvents, 1));

% Load State Config (Colors/Names)
try
    cfg = as_loadConfig([]);
    colors = cfg.colors;
    stateNames = cfg.names;
catch
    colors = lines(nStates);
    stateNames = arrayfun(@(x) sprintf('State %d', x), ...
        1:nStates, 'UniformOutput', false);
end

% Ensure colors is a cell array for processing consistency
if isnumeric(colors)
    colors = num2cell(colors, 2);
end

% =========================================================================
%  CALCULATION (PER BOUT)
%  ========================================================================

evtDur = evtTimes(:,2) - evtTimes(:,1);
tblState = cell(nStates, 1);
stateCounts = zeros(nStates, 1);

for iState = 1:nStates
    bouts = boutTimes{iState};
    nBouts = size(bouts, 1);

    % Skip table generation if state is empty
    if nBouts == 0
        continue;
    end

    % Map events to bouts for this state (Categorical Assignment)
    if ~isempty(bouts)
        inState = intervals(bouts).contains(peakTimes);
        stateIdx(inState) = stateNames{iState};
        stateCounts(iState) = sum(inState & accepted);
    end

    % Pre-allocate vectors for this state
    rate = nan(nBouts, 1);
    density = nan(nBouts, 1);
    boutDur = bouts(:, 2) - bouts(:, 1);

    % Iterate Bouts to calculate Rate/Density
    for iBout = 1:nBouts
        tStart = bouts(iBout, 1);
        tEnd   = bouts(iBout, 2);

        % Find accepted events in this specific bout
        idxBout = peakTimes >= tStart & peakTimes <= tEnd & accepted;

        count = sum(idxBout);
        durSum = sum(evtDur(idxBout));

        % Rate (Hz) = Count / Bout Duration
        rate(iBout) = count / boutDur(iBout);

        % Density (fraction) = Total Event Duration / Bout Duration
        density(iBout) = durSum / boutDur(iBout);
    end

    % Create Table for this state
    tblState{iState} = table(rate, density, boutDur, ...
        repmat(categorical(stateNames(iState)), nBouts, 1), ...
        bouts(:,1), bouts(:,2), ...
        'VariableNames', {'Rate', 'Density', 'Duration', 'State', 'Start', 'End'});
end

% Combine all states into one master table (empty when there are no bouts)
if isempty(tblState)
    evtStates = table();
else
    evtStates = vertcat(tblState{:});
end

% =========================================================================
%  OUTPUT & SAVING
%  ========================================================================

if flgSave
    % save the bout table under the modality-named variable so existing
    % readers (e.g. mcu_tblVivo expects 'evtStates') keep working
    outStruct.([name, 'States']) = evtStates;
    save(savefile, '-struct', 'outStruct', '-v7.3');
end

% =========================================================================
%  PLOTTING
%  ========================================================================
if flgPlot && ~isempty(evtStates)

    fh = figure('Name', [basename '_' name 'States'], 'NumberTitle', 'off');
    tiledlayout(2, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

    T = evtStates;
    plotColors = vertcat(colors{:});

    % Rate vs Time
    nexttile([1, 2]); hold on;
    gscatter(T.Start / 3600, T.Rate, T.State, plotColors, '.', 10, 'off');
    xlabel('Time (h)'); ylabel('Rate (Hz)'); title([lbl ' Rate']);
    axis tight;

    % Density vs Time
    nexttile([1, 2]); hold on;
    gscatter(T.Start / 3600, T.Density, T.State, plotColors, '.', 10, 'off');
    xlabel('Time (h)'); ylabel('Density (s/s)'); title([lbl ' Density']);
    axis tight;

    % Rate vs Duration (Check for short-bout bias)
    nexttile([1, 2]); hold on;
    gscatter(T.Duration, T.Rate, T.State, plotColors, '.', 10, 'off');
    xlabel('Bout Duration (s)'); ylabel('Rate (Hz)'); title('Rate vs Duration');
    set(gca, 'XScale', 'log'); axis tight;

    % Pie Chart (Total Counts)
    nexttile([1, 2]);
    if sum(stateCounts) > 0
        hPie = pie(stateCounts);

        % Apply Colors to Pie Chart
        % pie returns text/patch objects. Patches are at odd indices (1,3,5...)
        cnt = 1;
        for k = 1:2:length(hPie)
            if cnt <= length(colors)
                hPie(k).FaceColor = colors{cnt};
                cnt = cnt + 1;
            end
        end

        legend(stateNames, 'Location', 'bestoutside');
        title([lbl ' Count Distribution']);
    end

    sgtitle([basename ' - ' lbl ' States'], 'Interpreter', 'none');

    figDir = fullfile(basepath, 'graphics');
    if ~exist(figDir, 'dir'), mkdir(figDir); end
    saveas(fh, fullfile(figDir, [basename '_' name '_states.png']));
end

end