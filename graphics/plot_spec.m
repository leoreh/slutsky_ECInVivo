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
% 29 mar 22 LH      updates:
% 05 mar 24             changes mapping of x values according to tstamps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar)
addParameter(p, 'ch', 1, @isnumeric)
addParameter(p, 'logfreq', false, @islogical)
addParameter(p, 'saveFig', true, @(x) islogical(x) || ischar(x))
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

% params
freq = spec.freq;
nbins = length(spec.tstamps);
tspec = spec.tstamps / xtime;
nch = length(ch);

% file
[~, basename] = fileparts(basepath);

% open figure if no axis provided   
if isempty(axh)
    clear axh
    fh = figure;
    set(fh, 'WindowState', 'maximized');
    tlayout = [nch, 1];
    th = tiledlayout(nch, 1);
    th.TileSpacing = 'tight';
    th.Padding = 'none';
    title(th, basename, 'interpreter', 'none')
end

for ich = 1 : nch

    if exist('fh', 'var')
        axh(ich) = nexttile(th, ich, [1, 1]); cla; hold on
    end
    
    s = spec.s(:, :, ch(ich));

    % take a sample of the spectrogram to initialize the colormap (robust to
    % short inputs and to degenerate / non-finite ranges)
    nSamp = min(length(tspec), max(1, round(length(tspec) / 5)));
    sampleBins = randperm(length(tspec), nSamp);
    specSample = reshape(s(sampleBins, :), 1, []);
    cLimit = prctile(specSample, [2 94]);
    if numel(cLimit) < 2 || ~all(isfinite(cLimit)) || cLimit(1) >= cLimit(2)
        cLimit = [min(s(:)), max(s(:))];
        if ~all(isfinite(cLimit)) || cLimit(1) >= cLimit(2)
            cLimit = [0 1];
        end
    end
    
    imagesc(axh(ich), tspec, freq, s', cLimit);
    colormap(axh(ich), AccuSleep_colormap());
    axis(axh(ich), 'xy');
    ylabel(axh(ich), 'Frequency [Hz]');
    xlabel(axh(ich), 'Time [h]');
    if logfreq
        set(axh(ich), 'yscale', 'log');
        ylim(axh, [max([0.2, freq(1)]), freq(end)]);
    else
        set(axh(ich), 'yscale', 'linear');
    end

end

if length(axh) > 1
    linkaxes(axh, 'x')
end
axis(axh, 'tight')

% save (only when this function created its own figure; when an axis handle
% is supplied there is no figure to save and we must not call gca/savefig)
if saveFig && exist('fh', 'var')

    if ischar(saveFig)
        basepath = saveFig;
    end

    [~, basename] = fileparts(basepath);
    figpath = fullfile(basepath, 'graphics');
    if ~exist(figpath, 'dir'), mkdir(figpath); end
    figname = fullfile(figpath, sprintf('%s_spec', basename));
    savefig(fh, figname)
end

end

% EOF
