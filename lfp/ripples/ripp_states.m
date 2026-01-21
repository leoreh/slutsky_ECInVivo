function rippStates = ripp_states(rippTimes, peakTimes, boutTimes, varargin)
% RIPP_STATES Analyzes ripple events in relation to vigilance states.
%
% SUMMARY:
% This function calculates ripple rates and densities within different
% vigilance states.
%
% INPUT:
%   rippTimes           [N x 2] Ripple start/end times.
%   peakTimes            [N x 1] Ripple peak times.
%   boutTimes           Cell array of state times (e.g. ss.bouts.times).
%   basepath            (Optional) Path to session {pwd}.
%   flgPlot             (Optional) Plot results {true}.
%   flgSave             (Optional) Save 'rippStates' struct to disk {true}.
%   flgSaveFig          (Optional) Save figure {true}.
%
% OUTPUT:
%   rippStates          Structure with rate, density, and indices.
%   * Saves 'basename.rippStates.mat'.
%
% DEPENDENCIES: basepaths2vars, InIntervals, times2rate, as_loadConfig, setMatlabGraphics.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'rippTimes', @isnumeric);
addRequired(p, 'peakTimes', @isnumeric);
addRequired(p, 'boutTimes', @iscell);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgSave', true, @islogical);
addParameter(p, 'flgSaveFig', true, @islogical);
parse(p, rippTimes, peakTimes, boutTimes, varargin{:});

rippTimes = p.Results.rippTimes;
peakTimes = p.Results.peakTimes;
boutTimes = p.Results.boutTimes;
basepath = p.Results.basepath;
flgPlot = p.Results.flgPlot;
flgSave = p.Results.flgSave;
flgSaveFig = p.Results.flgSaveFig;

%% ========================================================================
%  INITIALIZATION
%  ========================================================================
cd(basepath);
[~, basename] = fileparts(basepath);
savefile = fullfile(basepath, [basename, '.rippStates.mat']);

nStates = length(boutTimes);
nRipples = length(peakTimes);

% Determine recording window from events if not passed (or just use global min/max)
% But for rate calculation usually we need recWin.
% However, prior code used `ripp.info.recWin`. Now we don't have it.
% We will use min/max of boutTimes or rippTimes as effective window?
% Or just assume [0, max(end)].
% Actually, simpler: define recWin based on max of everything.
maxTime = max([max(rippTimes(:)); 10]);
for i=1:nStates
    if ~isempty(boutTimes{i})
        maxTime = max(maxTime, max(boutTimes{i}(:)));
    end
end
recWin = [0 maxTime];

% Initialize Output Struct
rippStates = struct();
rippStates.rate = cell(nStates, 1);     % Rate [Hz] per bout
rippStates.density = cell(nStates, 1);  % Density [s/s] per bout
rippStates.idx = false(nRipples, nStates); % [N x S] Logical index
rippStates.binedges = cell(nStates, 1);
rippStates.tstamps = cell(nStates, 1);

%% ========================================================================
%  CALCULATION
%  ========================================================================
fprintf('Calculating ripple states for %s...\n', basename);

% Calculate Duration
dur = rippTimes(:,2) - rippTimes(:,1);

for iState = 1:nStates
    bouts = boutTimes{iState};

    if isempty(bouts)
        rippStates.rate{iState} = NaN;
        rippStates.density{iState} = NaN;
        continue;
    end

    % Ripple Rate (per bout)
    [rippStates.rate{iState}, rippStates.binedges{iState}, rippStates.tstamps{iState}] = ...
        times2rate(peakTimes, 'binsize', Inf, 'winCalc', bouts, 'c2r', true);

    % Ripple Density (per bout)
    nBouts = size(bouts, 1);
    boutDensity = zeros(nBouts, 1);
    for iBout = 1:nBouts
        inBout = (peakTimes >= bouts(iBout,1)) & (peakTimes <= bouts(iBout,2));

        % Density = Total Duration of Ripples / Bout Duration
        totalRippDur_s = sum(dur(inBout)); % Times are in seconds usually?
        % Wait, in wrapper 'ripp.dur' was converted to ms?
        % ripp_detect returns dur in ms. Wrapper calculates rippSamps from ripp.times (sec).
        % ripp_states takes rippTimes (sec).
        % ripp.dur in detect is converted to ms at end.
        % So here rippTimes is likely seconds.
        % Dur calculation: rippTimes(:,2)-rippTimes(:,1) is seconds.
        % So totalRippDur_s is sum(dur).

        boutDur_s = bouts(iBout, 2) - bouts(iBout, 1);

        if boutDur_s > 0
            boutDensity(iBout) = totalRippDur_s / boutDur_s;
        end
    end
    rippStates.density{iState} = boutDensity;

    % Index
    rippStates.idx(:, iState) = InIntervals(peakTimes, bouts);
end

%% ========================================================================
%  OUTPUT & SAVING
%  ========================================================================

if flgSave
    save(savefile, 'rippStates', '-v7.3');
    fprintf('Saved: %s\n', savefile);
end

%% ========================================================================
%  PLOTTING
%  ========================================================================
if flgPlot

    stateNames = arrayfun(@(x) sprintf('State %d',x), 1:nStates, 'UniformOutput', false);

    % Attempt to load colors config
    try
        cfg = as_loadConfig([]);
        colors = cfg.colors;
    catch
        colors = lines(nStates);
    end

    setMatlabGraphics(true);
    fh = figure('Name', [basename '_rippleStates'], 'NumberTitle', 'off');

    % Rate
    subplot(2, 4, [1, 2]);
    hold on;
    for iState = 1:nStates
        if ~isempty(rippStates.tstamps{iState})
            plot(rippStates.tstamps{iState} / 3600, ...
                rippStates.rate{iState}, '.', ...
                'Color', colors{iState}, ...
                'MarkerSize', 10);
        end
    end
    xlabel('Time (h)'); ylabel('Rate (Hz)'); title('Ripple Rate');
    axis tight;

    % Density
    subplot(2, 4, [5, 6]);
    hold on;
    for iState = 1:nStates
        if ~isempty(rippStates.tstamps{iState})
            plot(rippStates.tstamps{iState} / 3600, ...
                rippStates.density{iState}, '.', ...
                'Color', colors{iState}, ...
                'MarkerSize', 10);
        end
    end
    xlabel('Time (h)'); ylabel('Density (s/s)'); title('Ripple Density');
    axis tight;

    % Pie
    subplot(2, 4, [3, 4, 7, 8]);
    counts = sum(rippStates.idx, 1);
    if sum(counts) > 0
        pie(counts, ones(size(counts)));
        legend(stateNames, 'Location', 'bestoutside');
        title('Ripple State Distribution');
    end

    sgtitle([basename ' - Ripple States'], 'Interpreter', 'none');

    if flgSaveFig
        figDir = fullfile(basepath, 'graphics');
        if ~exist(figDir, 'dir')
            mkdir(figDir);
        end
        saveas(fh, fullfile(figDir, [basename '_ripp_states.png']));
    end
end

end