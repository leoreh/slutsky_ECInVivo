function [stateIdx, rippStates] = ripp_states(rippTimes, peakTimes, boutTimes, varargin)
% RIPP_STATES Classifies ripples by vigilance state and computes bout statistics.
%
%   [stateIdx, rippStates] = RIPP_STATES(rippTimes, peakTimes, boutTimes, varargin)
%
%   SUMMARY:
%       1. Assigns a vigilance state to each ripple based on its peak time.
%       2. Calculates Ripple Rate (Hz) and Density (duration/duration) for
%          every individual sleep bout.
%       3. Aggregates this into a summary table.
%
%   INPUTS:
%       rippTimes   - (Mat)  [N x 2] Ripple start/end times [s].
%       peakTimes   - (Vec)  [N x 1] Ripple peak times [s].
%       boutTimes   - (Cell) {N_states x 1} Cell array of [Start End] matrices.
%                            Typically {Wake, NREM, REM, ...}.
%       varargin    - Parameter/Value pairs:
%           'basepath' - (Char) Save location. (Default: pwd).
%           'flgPlot'  - (Log)  Generate summary figures? (Default: true).
%           'flgSave'  - (Log)  Save .rippStates.mat? (Default: true).
%
%   OUTPUTS:
%       stateIdx    - (Cat)   [N_ripples x 1] Categorical array of states.
%       rippStates  - (Table) Summary table with rows per bout:
%                             [Rate, Density, Duration, State, Start, End].
%
%   DEPENDENCIES:
%       basepaths2vars, InIntervals, as_loadConfig.
%
%   HISTORY:
%       Updated: 23 Jan 2026

% =========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'rippTimes', @isnumeric);
addRequired(p, 'peakTimes', @isnumeric);
addRequired(p, 'boutTimes', @iscell);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgSave', true, @islogical);
parse(p, rippTimes, peakTimes, boutTimes, varargin{:});

rippTimes = p.Results.rippTimes;
peakTimes = p.Results.peakTimes;
boutTimes = p.Results.boutTimes;
basepath = p.Results.basepath;
flgPlot = p.Results.flgPlot;
flgSave = p.Results.flgSave;

% =========================================================================
%  PREP
%  ========================================================================
cd(basepath);
[~, basename] = fileparts(basepath);
savefile = fullfile(basepath, [basename, '.rippStates.mat']);

nStates = length(boutTimes);
nRipples = length(peakTimes);

% Initialize Categorical Index
stateIdx = categorical(nan(nRipples, 1));

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

rippDur = rippTimes(:,2) - rippTimes(:,1);
tblState = cell(nStates, 1);
stateCounts = zeros(nStates, 1);

for iState = 1:nStates
    bouts = boutTimes{iState};
    nBouts = size(bouts, 1);

    % Skip table generation if state is empty
    if nBouts == 0
        continue;
    end

    % Map ripples to bouts for this state (Categorical Assignment)
    if ~isempty(bouts)
        inState = InIntervals(peakTimes, bouts);
        stateIdx(inState) = stateNames{iState};
        stateCounts(iState) = sum(inState);
    end

    % Pre-allocate vectors for this state
    rate = nan(nBouts, 1);
    density = nan(nBouts, 1);
    boutDur = bouts(:, 2) - bouts(:, 1);

    % Iterate Bouts to calculate Rate/Density
    for iBout = 1:nBouts
        tStart = bouts(iBout, 1);
        tEnd   = bouts(iBout, 2);

        % Find ripples in this specific bout
        idxBout = peakTimes >= tStart & peakTimes <= tEnd;

        count = sum(idxBout);
        durSum = sum(rippDur(idxBout));

        % Rate (Hz) = Count / Bout Duration
        rate(iBout) = count / boutDur(iBout);

        % Density (fraction) = Total Ripple Duration / Bout Duration
        density(iBout) = durSum / boutDur(iBout);
    end

    % Create Table for this state
    tblState{iState} = table(rate, density, boutDur, ...
        repmat(categorical(stateNames(iState)), nBouts, 1), ...
        bouts(:,1), bouts(:,2), ...
        'VariableNames', {'Rate', 'Density', 'Duration', 'State', 'Start', 'End'});
end

% Combine all states into one master table
rippStates = vertcat(tblState{:});

% =========================================================================
%  OUTPUT & SAVING
%  ========================================================================

if flgSave
    save(savefile, 'rippStates', '-v7.3');
    if flgSave
        save(savefile, 'rippStates', '-v7.3');
    end
end

% =========================================================================
%  PLOTTING
%  ========================================================================
if flgPlot

    fh = figure('Name', [basename '_rippleStates'], 'NumberTitle', 'off');
    tiledlayout(2, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

    T = rippStates;
    plotColors = vertcat(colors{:});

    % Rate vs Time
    nexttile([1, 2]); hold on;
    gscatter(T.Start / 3600, T.Rate, T.State, plotColors, '.', 10, 'off');
    xlabel('Time (h)'); ylabel('Rate (Hz)'); title('Ripple Rate');
    axis tight;

    % Density vs Time
    nexttile([1, 2]); hold on;
    gscatter(T.Start / 3600, T.Density, T.State, plotColors, '.', 10, 'off');
    xlabel('Time (h)'); ylabel('Density (s/s)'); title('Ripple Density');
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
        title('Ripple Count Distribution');
    end

    sgtitle([basename ' - Ripple States'], 'Interpreter', 'none');

    figDir = fullfile(basepath, 'graphics');
    if ~exist(figDir, 'dir'), mkdir(figDir); end
    saveas(fh, fullfile(figDir, [basename '_ripp_states.png']));
end

end