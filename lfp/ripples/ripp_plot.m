function ripp_plot(ripp, varargin)
% RIPP_PLOT Generates summary plots for ripple detection results.
%
% SUMMARY:
% This function takes the output structure from ripp_detect.m and creates
% a figure summarizing key ripple properties, including:
%   - Overall ripple rate over time.
%   - Average ripple waveform (filtered LFP).
%   - Peri-event maps of frequency and amplitude.
%   - Ripple autocorrelogram (ACG).
%   - Distributions of peak frequency, peak amplitude, and duration.
%
% INPUT:
%   ripp          	Structure containing ripple detection results from
%                   ripp_detect.m. Must include fields like ripp.rate, 
%                   ripp.maps, ripp.acg, ripp.peakFreq, ripp.peakAmp, ripp.dur.
%   basepath        (Optional) Path to recording session directory {pwd}.
%   flgSaveFig      (Optional) Logical flag to save the figure {true}.
%
% OUTPUT:
%   None. Generates and optionally saves a figure.
%
% DEPENDENCIES:
%   setMatlabGraphics (custom), plot_stdShade (custom), PlotColorMap (custom), 
%   plot_ccg (custom), export_fig (external toolbox)
%
% HISTORY:
% Aug 2024 LH - Refactored from plot_ripples.m to align with ripp_detect style.
% 11 Jan 23 LH - Original version (plot_ripples).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addRequired(p, 'ripp', @isstruct); % Require the ripple structure
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSaveFig', true, @islogical);

parse(p, ripp, varargin{:}); % Pass 'ripp' directly
ripp        = p.Results.ripp; % Retrieve ripp from parsed results
basepath    = p.Results.basepath;
flgSaveFig  = p.Results.flgSaveFig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARATIONS & PARAMETER EXTRACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Files
cd(basepath);
[~, basename] = fileparts(basepath);
rippfile = fullfile(basepath, [basename, '.ripp.mat']); % Define for consistency, though ripp is input

% Load ripp structure if empty (optional, useful for standalone execution)
if isempty(ripp)
    if exist(rippfile, 'file')
        load(rippfile, 'ripp')
    else
        error('ripp file missing')
    end
end


% Parameters for plotting
nRipples = size(ripp.times, 1);
mapDur = ripp.maps.durWin; % Use map window from ripp struct [s]
nbinsMap = size(ripp.maps.freq, 2);
xMap = linspace(mapDur(1), mapDur(2), nbinsMap); % Time axis for maps [s]
histBins = 100; % Number of bins for histograms
plotNumRipples = min([100, nRipples]); % Number of ripples for maps
ripp_idx = randperm(nRipples, plotNumRipples); % Random subset for visualization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHICS GENERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMatlabGraphics(true); % Apply consistent graphic settings
fh = figure('Name', [basename '_ripples'], 'NumberTitle', 'off');

% --- Rate Plot ---
sb1 = subplot(3, 3, [1, 2]);
plot(ripp.rate.timestamps / 3600, ripp.rate.rate, 'k', 'LineWidth', 1); % Plot rate in hours
xlabel('Time [h]');
ylabel('Ripple Rate [Hz]');
title('Overall Ripple Rate');
axis tight;

% --- Average Ripple Waveform ---
sb2 = subplot(3, 3, 3);
plot_stdShade('axh', sb2, 'dataMat', ripp.maps.filt', 'alpha', 0.3,...
    'clr', [0 0 0], 'xVal', xMap); % Mean +/- std of filtered LFP
xlabel('Time relative to peak [s]');
ylabel('Filtered LFP (norm.)'); % Assuming maps are normalized or z-scored
title('Average Ripple Shape');
xlim(mapDur);

% --- Frequency Map ---
sb3 = subplot(3, 3, 4);
PlotColorMap(ripp.maps.freq(ripp_idx, :), 'x', xMap, 'bar','on', ...
    'cutoffs', ripp.info.passband); % Passband from ripp.info
ylabel(sprintf('Ripple # (rand %d)', plotNumRipples));
xlabel('Time relative to peak [s]');
title('Instantaneous Frequency');
xlim(mapDur);

% --- Amplitude Map ---
sb4 = subplot(3, 3, 5);
PlotColorMap(ripp.maps.amp(ripp_idx, :), 'x', xMap, 'bar','on', 'cutoffs', []);
ylabel(sprintf('Ripple #'));
xlabel('Time relative to peak [s]');
title('Amplitude Envelope');
xlim(mapDur);

% --- Autocorrelogram (ACG) ---
sb5 = subplot(3, 3, 6);
plot_ccg(ripp.acg.data, ripp.acg.t * 1000); % Convert ACG time to ms
xlabel('Time Lag [ms]');
ylabel('Counts');
title('Ripple ACG');

% --- Peak Frequency Distribution ---
sb6 = subplot(3, 3, 7);
h = histogram(ripp.peakFreq, histBins, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1.5, 'EdgeColor', 'k');
xlabel('Peak Frequency [Hz]');
ylabel('Probability');
title('Peak Frequency Dist.');
box off;

% --- Peak Amplitude Distribution ---
sb7 = subplot(3, 3, 8);
h = histogram(ripp.peakAmp, histBins, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1.5, 'EdgeColor', 'k');
% Consider log scale if distribution is skewed, but start linear
% set(gca, 'xscale', 'log'); 
xlabel('Peak Amplitude (norm.)'); % Adjust label if not normalized
ylabel('Probability');
title('Peak Amplitude Dist.');
box off;

% --- Duration Distribution ---
sb8 = subplot(3, 3, 9);
h = histogram(ripp.dur, histBins, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1.5, 'EdgeColor', 'k');
set(gca, 'xscale', 'log'); % Duration often benefits from log scale
xlabel('Duration [ms]');
ylabel('Probability');
title('Duration Dist.');
box off;
xticks([10, 20, 50, 100, 200, 300]); % Example sensible ticks for log scale

% --- Overall Title ---
sgtitle([basename ' - Ripple Summary'], 'Interpreter', 'none');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flgSaveFig
    figpath = fullfile(basepath, 'graphics');
    if ~exist(figpath, 'dir')
        mkdir(figpath);
    end
    figname = fullfile(figpath, sprintf('%s_ripp_plot', basename)); % Use new function name
    saveas(fh, [figname '.png']); 
end

end

% EOF