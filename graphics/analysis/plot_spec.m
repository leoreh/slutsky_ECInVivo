function plot_spec(spec, logfreq, saveFig)

% plots a spectrogram on the current figure using the struct created in
% calc_spec.
%
% INPUT
%   spec        struct. see calc_spec
%   logfreq     logical. plot y-axis on log (true) or linear {false} scale
%   saveFig     logical / char. if char will be treated as the filepath for
%               saving the figure
%
% 29 mar 22 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% laze fix for multichannel
if ndims(spec.s) == 3
    spec.s = squeeze(spec.s(:, :, 1));
end

% params
winstep = spec.info.winstep;
freq = spec.freq;

% time axis in hours
nbins = length(spec.tstamps);
tspec = ((1 : nbins) * winstep - winstep / 2) / 3600;

% take a sample of the spectrogram to help initialize the colormap
sampleBins = randperm(nbins, round(nbins / 10));

% plot
if ~logfreq
    
    specSample = reshape(spec.s(sampleBins, :), 1, length(sampleBins) * length(freq));
    cLimit = prctile(specSample, [6 98]);

    imagesc(tspec, freq, spec.s', cLimit);
else
    pband = 10 * log10(abs(spec.s));
    pband = bz_NormToRange(pband, [0 1]);
    surf(tspec, freq, pband', 'EdgeColor', 'none');
    set(gca, 'yscale', 'log')
    view(0, 90);
    ylim([max([0.5, freq(1)]), freq(end)])

    specSample = reshape(pband(sampleBins, :), 1, length(sampleBins) * length(freq));
    cLimit = prctile(specSample, [35 99.9]);
    clim(cLimit)
end
colormap(AccuSleep_colormap());
axis('xy')
ylabel('Frequency [Hz]')
xlabel('Time [h]')

% save
if saveFig
    
    if ischar(saveFig)
        basepath = saveFig;
    else
        basepath = pwd;
    end
    
    [~, basename] = fileparts(basepath);
    figpath = fullfile(basepath, 'graphics');
    mkdir(figpath)
    figname = fullfile(figpath, sprintf('%s_spec', basename));
    export_fig(figname, '-tif', '-transparent', '-r300')
end

end

% EOF

