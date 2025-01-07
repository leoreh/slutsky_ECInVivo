function plot_ripples(ripp, varargin)

% plots spiking activity in realtion to ripples. see getRipples and
% getRippleSpks.m
%
% INPUT:
%   ripp          	struct
%   basepath        path to recording {pwd}
%   saveFig         logical
%
% 11 jan 23 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'saveFig', true, @islogical);

parse(p, varargin{:})
basepath    = p.Results.basepath;
saveFig     = p.Results.saveFig;

% files
cd(basepath)
[~, basename] = fileparts(basepath);
rippfile = fullfile(basepath, [basename, '.ripp.mat']);

% load
if isempty(ripp)
    if exist(rippfile, 'file')
        load(rippfile, 'ripp')
    else
        error('ripp file missing')
    end
end

% params
nbouts = size(ripp.bouts, 1);
durPlot = [-50 50] / 1000;
x = durPlot(1) : diff(durPlot) / nbouts : durPlot(2);
histBins = 200;
nbinsMap = size(ripp.maps.freq, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMatlabGraphics(true)
fh = figure;

% rate
sb1 = subplot(3, 3, [1, 2]);
plot(ripp.rate.tstamps / 60 / 60, ripp.rate.rate, 'k')
xlabel('Time [h]')
ylabel('Ripple Rate [1/min]')

% examples of ripples
ripp_idx = randperm(nbouts, min([100, nbouts]));

% mean +/- std of ripples (filtered) superimposed
sb2 = subplot(3, 3, 3);
xval = ((1 : nbinsMap)' - ceil(nbinsMap / 2)) / nbinsMap * diff(durPlot);
plot_stdShade(sb2, ripp.maps.ripp, 0.3, 'k', xval, [])
xlabel('Time [s]')

% frequency map
sb3 = subplot(3, 3, 4);
PlotColorMap(ripp.maps.freq(ripp_idx, :), 1, 'bar','on',...
    'cutoffs', ripp.info.passband, 'x', x);
ylabel('Ripple No.')
title('Frequency');

% amplitude map
sb4 = subplot(3, 3, 5);
PlotColorMap(ripp.maps.amp(ripp_idx, :), 1, 'bar','on', 'x', x, 'cutoffs', []);
ylabel('Ripple No.')
title('Amplitude');

% ACG
sb5 = subplot(3, 3, 6);
plot_ccg(ripp.acg.data, ripp.acg.t);
xlabel('Time [ms]')
ylabel('Rate')
title('Auto correlation');

% distribution of peak frequency
sb6 = subplot(3, 3, 7);
h = histogram(ripp.peakFreq, histBins, 'Normalization', 'probability');
h.FaceColor = 'k';
h.EdgeColor = 'none';
xlabel('Peak Frequency [Hz]')
ylabel('Probability')

% distribution of peak amplitude
sb7 = subplot(3, 3, 8);
h = histogram(ripp.peakAmp, histBins, 'Normalization', 'probability');
h.FaceColor = 'k';
h.EdgeColor = 'none';
xlabel('Peak Amp')
ylabel('Probability')
set(gca, 'xscale', 'log')

% distribution of ripple duration
sb8 = subplot(3, 3, 9);
h = histogram(ripp.dur * 1000, histBins, 'Normalization', 'probability');
h.FaceColor = 'k';
h.EdgeColor = 'none';
xlabel('Ripple Duration [log(ms)]')
ylabel('Probability')
set(gca, 'xscale', 'log')
xticks([1, 10, 20, 50, 100, 200, 300])

sgtitle(basename)

% save figure
if saveFig
    figpath = fullfile(basepath, 'graphics');
    mkdir(figpath)
    figname = fullfile(figpath, sprintf('%s_ripples', basename));
    export_fig(figname, '-tif', '-transparent', '-r300')
end

end

% EOF