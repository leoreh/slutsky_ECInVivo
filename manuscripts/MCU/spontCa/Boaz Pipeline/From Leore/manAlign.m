    function manShift = manAlign(varargin)

% plots two signals for visual alignment. The 1st signal is spktimes and
% is plotted as a raster plot. The 2nd is a time series. 
%
% INPUT:
%   spktimes    numeric vector of spktimes [s].
%   sig         numeric vector of time series. 
%   fs          sampling frequency of sig
%
% DEPENDENCIES
%   plotSpikeRaster
%
%  15 may 23 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'spktimes', [], @isnumeric);
addOptional(p, 'sig', [], @isnumeric);
addOptional(p, 'fs', [], @isnumeric);

parse(p, varargin{:})
spktimes    = p.Results.spktimes;
sig         = p.Results.sig;
fs          = p.Results.fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% timestamps for sig
tstamps = [1 : length(sig)] / fs;
recDur = tstamps(end);

% viewing parameters
zoomFactor = 0.1;                   % zoom in factor for second plot
xStart = round(recDur / 2);         % time to start second plot [s]
xRange = (floor(length(sig) * zoomFactor)) / fs;

% raster plot params
LineFormat.Color = [0.8 0.1 0.1];

% shift var
global manShift

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a figure
fh = figure();

% 1st plot of entire signal
axh1 = subplot(4, 1, 1);
plotSpikeRaster({spktimes'}, 'PlotType', 'vertline', 'LineFormat', LineFormat);
axh2 = subplot(4, 1, 2);
pf1 = plot(tstamps, sig);
linkaxes([axh1, axh2], 'x')

% 2nd plot for zoom in view of center segment
axh3 = subplot(4, 1, 3);
plotSpikeRaster({spktimes'}, 'PlotType', 'vertline', 'LineFormat', LineFormat);
axh4 = subplot(4, 1, 4);
pf2 = plot(tstamps, sig);
yLimit = ylim;
set(gca, 'YLimMode', 'manual');
ylim(yLimit)
linkaxes([axh3, axh4], 'x')
xlim([xStart, xStart + xRange])

% Create slider for adjusting the shift
slider = uicontrol('Style', 'slider', 'Units', 'normalized',...
    'Position', [0.1, 0.05, 0.8, 0.03], 'Value', 0, 'Min', -0.5, 'Max', 0.5);

% define slider resolution (minor and major step size)
stepSize = 0.001; 
slider.SliderStep = [stepSize, stepSize * 10]; 

% add callback to update the plot when the slider is adjusted
addlistener(slider, 'Value', 'PostSet', @(~, ~) updatePlot());

% update the plot based on the slider value
    function updatePlot()
        % Get the current value of the slider
        manShift = slider.Value;

        % Shift the second signal by the specified amount
        sig_shifted = circshift(sig, round(manShift * length(tstamps)));

        % Update the plot of the shifted signal
        set(pf1, 'YData', sig_shifted);
        set(pf2, 'YData', sig_shifted);

        % Refresh the figure
        drawnow;
    end

uiwait(fh)
manShift = manShift * recDur;

end