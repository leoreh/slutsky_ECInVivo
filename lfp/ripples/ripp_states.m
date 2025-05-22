function ripp = ripp_states(ripp, varargin)
% RIPP_STATES Analyzes ripple events in relation to vigilance states.
%
% SUMMARY:
% This function takes ripple event data and vigilance state information 
% (from sleep_states.mat via ss_classify.m) to characterize ripple 
% occurrence within different states. It calculates:
%   1. Ripple rate within each defined vigilance state.
%   2. Ripple density (total duration per bout) within each state.
%   3. An index indicating which ripples occurred during each state.
%   4. Stores these results in the input 'ripp' structure.
%   5. Optionally plots the overall ripple rate overlaid with state-specific 
%      rates/densities and a pie chart showing the distribution of ripples across states.
%   6. Optionally saves the updated 'ripp' structure.
%
% NOTE ON CALCULATIONS:
% For state-specific rates, we calculate one rate value per state bout by
% setting binsize to Inf. This approach ensures statistical independence
% between data points and avoids temporal autocorrelation within bouts.
%
% For state-specific densities, we calculate one density value per bout by
% summing the total duration of ripples that occur within that bout and
% dividing by the bout duration. This provides a measure of how much time
% is spent in ripples during each bout.
%
% INPUT:
%   ripp                Structure containing ripple detection results. 
%                       Must include fields: 
%                         ripp.peakTime: [N x 1] ripple peak times.
%                         ripp.times: [N x 2] ripple start/end times.
%                         ripp.dur: [N x 1] ripple durations in ms.
%                         ripp.info.fs: LFP sampling frequency.
%                         ripp.info.recWin: [1 x 2] recording window used.
%                         ripp.info.binsizeRate: Binsize for rate calc [s].
%                         ripp.rate.rate: Overall ripple rate vector.
%                         ripp.rate.tstamps: Timestamps for overall rate.
%   basepath            (Optional) Path to recording session directory {pwd}.
%   flgGraphics         (Optional) Logical flag to plot results {true}.
%   flgSaveVar          (Optional) Logical flag to save/update the ripp
%                       structure in a .ripp.mat file {true}.
%   flgSaveFig          (Optional) Logical flag to save the plot {true}.
%
% OUTPUT:
%   ripp                Updated structure containing ripple analysis results,
%                       with the following fields added or updated under
%                       ripp.states:
%                       ripp.states.stateNames: Cell array of state names.
%                       ripp.states.rate: {nStates x 1} cell array, each cell 
%                                         containing the ripple rate vector 
%                                         within that state [Hz].
%                       ripp.states.density: {nStates x 1} cell array, each cell 
%                                           containing the ripple density vector 
%                                           within that state [s/s].
%                       ripp.states.binedges: {nStates x 1} cell array, bin 
%                                             edges for state calc [s].
%                       ripp.states.tstamps: {nStates x 1} cell array, 
%                                            timestamps for state calc [s].
%                       ripp.states.idx: [nEvents x nStates] logical matrix 
%                                        indicating if a ripple (row) occurred 
%                                        during a state (column).
%
% DEPENDENCIES:
%   basepaths2vars (custom)
%   InIntervals (FMAToolbox) 
%   times2rate (custom)
%   as_loadConfig (custom, for plotting colors)
%   setMatlabGraphics (custom, for plotting)
%   export_fig (external toolbox, for saving figures)
%   Requires FMAToolbox in path for InIntervals.
%
% HISTORY:
% Aug 2024 LH - Major refactor to align with ripp_detect/ripp_spks style.
%               Renamed from rippleStates to ripp_states.
% 12 Jan 23 LH - Original version (rippleStates).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addRequired(p, 'ripp', @isstruct); % Require the ripple structure as input
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgGraphics', true, @islogical);
addParameter(p, 'flgSaveVar', true, @islogical);
addParameter(p, 'flgSaveFig', true, @islogical);

parse(p, ripp, varargin{:}); % Pass 'ripp' directly to parse
ripp                = p.Results.ripp; % Retrieve ripp from parsed results
basepath            = p.Results.basepath;
flgGraphics         = p.Results.flgGraphics;
flgSaveVar          = p.Results.flgSaveVar;
flgSaveFig          = p.Results.flgSaveFig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARATIONS & DATA LOADING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Files
cd(basepath);
[~, basename] = fileparts(basepath);
rippfile = fullfile(basepath, [basename, '.ripp.mat']);

% Load required session variables (sleep states)
vars = ["sleep_states"]; 
v = basepaths2vars('basepaths', {basepath}, 'vars', vars);

% Check if sleep states loaded
if ~isfield(v, 'ss') || ~isfield(v.ss.bouts, 'times') || ~isfield(v.ss, 'info') || ~isfield(v.ss.info, 'names')
    error('Failed to load required sleep state variables (ss.bouts.times, ss.info.names). Run ss_classify first.');
end
ss = v.ss;

% Parameters
peakTime    = ripp.peakTime;
recWin      = ripp.info.recWin;
nRipples    = length(peakTime);
nStates     = length(ss.bouts.times); % Number of states defined in ss

% Initialize output structure fields
ripp.states = []; % Clear any previous state analysis
ripp.states.stateNames = ss.info.names;
ripp.states.rate = cell(nStates, 1);
ripp.states.density = cell(nStates, 1);  % Add density field
ripp.states.binedges = cell(nStates, 1);
ripp.states.tstamps = cell(nStates, 1);
ripp.states.idx = false(nRipples, nStates); % Initialize logical index matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE RIPPLE PROPERTIES PER STATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iState = 1:nStates

    stateBouts = ss.bouts.times{iState};

    % Check if state bouts exist and overlap with the recording window
    if ~isempty(stateBouts)
        % Limit state bouts to the recording window analyzed for ripples
        validBoutIdx = InIntervals(mean(stateBouts, 2), recWin); % Check if bout center is within recWin
        stateBouts = stateBouts(validBoutIdx, :); 
        
        % Ensure bouts are strictly within recWin after selection (clip edges)
        stateBouts(:, 1) = max(stateBouts(:, 1), recWin(1));
        stateBouts(:, 2) = min(stateBouts(:, 2), recWin(2));
        stateBouts(diff(stateBouts, 1, 2) <= 0, :) = []; % Remove zero/negative duration bouts after clipping
    end

    % Proceed only if valid bouts remain within the recording window
    if ~isempty(stateBouts)
        % Calculate ripple rate restricted to these state bouts
        [ripp.states.rate{iState}, ripp.states.binedges{iState}, ...
         ripp.states.tstamps{iState}] = ...
            times2rate(peakTime, 'binsize', Inf, ...
                       'winCalc', stateBouts, 'c2r', true);

        % Calculate ripple density for each bout
        nBouts = size(stateBouts, 1);
        boutDensity = zeros(nBouts, 1);
        
        for iBout = 1:nBouts
            % Find ripples that occur during this bout
            boutRipples = InIntervals(peakTime, stateBouts(iBout,:));
            
            % Sum durations of ripples in this bout
            boutDensity(iBout) = sum(ripp.dur(boutRipples)) / (stateBouts(iBout,2) - stateBouts(iBout,1)) * 100 / 1000;
        end
        
        % Store density results
        ripp.states.density{iState} = boutDensity;

        % Identify ripples occurring within these state bouts
        ripp.states.idx(:, iState) = InIntervals(peakTime, stateBouts);

    else
        % If no valid bouts for this state, store NaNs/empty
        ripp.states.rate{iState} = NaN;
        ripp.states.density{iState} = NaN;
        ripp.states.binedges{iState} = NaN;
        ripp.states.tstamps{iState} = NaN;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE RESULTS & PLOT GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save updated ripp structure
if flgSaveVar
    save(rippfile, 'ripp', '-v7.3');
end

% Plotting
if flgGraphics
    
    % Use custom config for colors 
    cfg = as_loadConfig([]);
    plotColors = cfg.colors;
    sstates = 1 : length(ripp.states.rate);

    setMatlabGraphics(true); % Apply consistent graphic settings
    fh = figure('Name', [basename '_rippleStates'], 'NumberTitle', 'off');

    % Subplot 1: Ripple Rate over Time with State Overlays
    sb1 = subplot(2, 4, [1, 2]);
    plot(ripp.rate.timestamps / 3600, ripp.rate.rate, 'k', 'LineWidth', 1.5); % Plot overall rate in hours
    hold on;
    for iState = 1:nStates
        if ~isempty(ripp.states.tstamps{iState}) && ~all(isnan(ripp.states.rate{iState}))
            ph = plot(ripp.states.tstamps{iState} / 3600, ripp.states.rate{iState}, ...
                      '.', 'MarkerSize', 12, 'Color', plotColors{iState});
        end
    end
    xlabel('Time [h]');
    ylabel('Ripple Rate [Hz]');
    title('Ripple Rate');
    legend(['Overall'; ripp.states.stateNames], 'Location', 'best');
    hold off;
    axis tight;

    % Subplot 2: Ripple Density over Time with State Overlays
    sb2 = subplot(2, 4, [5, 6]);
    hold on;
    for iState = 1:nStates
        if ~isempty(ripp.states.tstamps{iState}) && ~all(isnan(ripp.states.density{iState}))
            ph = plot(ripp.states.tstamps{iState} / 3600, ripp.states.density{iState}, ...
                      '.', 'MarkerSize', 12, 'Color', plotColors{iState});
        end
    end
    xlabel('Time [h]');
    ylabel('Ripple Density [s/s]');
    title('Ripple Density');
    legend(ripp.states.stateNames, 'Location', 'best');
    hold off;
    axis tight;

    % Subplot 3: Pie Chart of Ripple Distribution across States
    sb3 = subplot(2, 4, [3, 4, 7, 8]);
    
    prct_states = sum(ripp.states.idx, 1);
    pie(prct_states, ones(1, length(prct_states)))
    hold on
    ph = findobj(sb3, 'Type', 'Patch');
    set(ph, {'FaceColor'}, flipud(cfg.colors(sstates)))
    legend({ripp.states.stateNames{sstates}}, 'FontSize', 10,...
        'Units', 'normalized', 'Location', 'best',...
        'NumColumns', 2);
      
    sgtitle([basename ' - Ripple Analysis by State'], 'Interpreter', 'none');

    % Save Figure
    if flgSaveFig
        figpath = fullfile(basepath, 'graphics');
        if ~exist(figpath, 'dir')
            mkdir(figpath);
        end
        figname = fullfile(figpath, sprintf('%s_ripp_states', basename));
        saveas(fh, [figname '.png']);
    end
end

end

% EOF