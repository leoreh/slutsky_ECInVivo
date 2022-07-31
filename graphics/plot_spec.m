function plot_spec(spec, varargin)

% plots a spectrogram on the current figure using the struct created in
% calc_spec. logfreq transforms the y-axis labels assuming that the spectrogram
% was sampled in log10 base, e.g. freq = logspace(log10(0.5), 2, 200). this
% is different then using surf (commented below)
%
% INPUT
%   spec        struct. see calc_spec
%   basepath    char. fullpath to recording folder {pwd}
%   ch          numeric. channel to plot from within spec
%   logfreq     logical. plot y-axis on log (true) or linear {false} scale
%   saveFig     logical / char. if char will be treated as the filepath for
%               saving the figure
%   xtime       numeric. factor by which to divide tstamps for x-axis.
%               e.g., 3600 will plot spec in hours
%   axh         axis handle for plot. if empty will create new figure
%
% 29 mar 22 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar)
addParameter(p, 'ch', 1, @isnumeric)
addParameter(p, 'logfreq', false, @islogical)
addParameter(p, 'saveFig', true, @islogical)
addParameter(p, 'xtime', 3600, @isnumeric)
addParameter(p, 'axh', [])

parse(p, varargin{:})
basepath        = p.Results.basepath;
ch              = p.Results.ch;
logfreq         = p.Results.logfreq;
saveFig         = p.Results.saveFig;
xtime           = p.Results.xtime;
axh             = p.Results.axh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lazy fix for multichannel spec
if ndims(spec.s) == 3
    spec.s = spec.s(:, :, ch);
end

% params
winstep = spec.info.winstep;
freq = spec.freq;

% time axis
nbins = length(spec.tstamps);
tspec = ((1 : nbins) * winstep - winstep / 2) / xtime;

% take a sample of the spectrogram to help initialize the colormap
sampleBins = randperm(nbins, round(nbins / 5));
specSample = reshape(spec.s(sampleBins, :), 1, length(sampleBins) * length(freq));
cLimit = prctile(specSample, [10 90]);

% plot   
if isempty(axh)
    fh = figure;
    axh = subplot(1, 1, 1);
end
imagesc(axh, tspec, freq, spec.s', cLimit);
colormap(AccuSleep_colormap());
axis('xy')
ylabel('Frequency [Hz]')
xlabel('Time [h]')
if logfreq
    set(gca, 'yscale', 'log')
    ylim([max([0.2, freq(1)]), freq(end)])
end

% save
if saveFig
    
    if ischar(saveFig)
        basepath = saveFig;
    end
    
    [~, basename] = fileparts(basepath);
    figpath = fullfile(basepath, 'graphics');
    mkdir(figpath)
    figname = fullfile(figpath, sprintf('%s_spec', basename));
    export_fig(figname, '-tif', '-transparent', '-r300')
end

end

% EOF

% if ~logfreq
%     
%     specSample = reshape(spec.s(sampleBins, :), 1, length(sampleBins) * length(freq));
%     cLimit = prctile(specSample, [6 98]);
% 
%     imagesc(tspec, freq, spec.s', cLimit);
% else
%     pband = 10 * log10(abs(spec.s));
%     pband = bz_NormToRange(pband, [0 1]);
%     surf(tspec, freq, pband', 'EdgeColor', 'none');
%     set(gca, 'yscale', 'log')
%     view(0, 90);
%     ylim([max([0.5, freq(1)]), freq(end)])
% 
%     specSample = reshape(pband(sampleBins, :), 1, length(sampleBins) * length(freq));
%     cLimit = prctile(specSample, [35 99.9]);
%     clim(cLimit)
% end
