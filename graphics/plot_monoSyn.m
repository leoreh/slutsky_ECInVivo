function plot_monoSyn(varargin)

% plot mono synaptic connection
% 
% INPUT:
%   basepath        path to recording {pwd}
%   spktimes        1 x 2 cell spike times [s]
%   units           1 x 2 numeric of unit idx
%   clr             2 element char. color for each cell
%   wv              2 x n numeric. mean waveform of units.
%   wv_std          2 x n numeric. std of units waveforms.
%   fs              numeric. sampling frequency 
%   saveFig         logical {true}
%
% DEPENDENCIES:
%   CCG
%
% 10 jan 22 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'spktimes', {}, @iscell);
addOptional(p, 'units', [], @isnumeric);
addOptional(p, 'clr', 'kk', @ischar);
addOptional(p, 'wv', [], @isnumeric);
addOptional(p, 'wv_std', [], @isnumeric);
addOptional(p, 'fs', 10000, @isnumeric);
addOptional(p, 'saveFig', true, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
spktimes        = p.Results.spktimes;
units           = p.Results.units;
clr             = p.Results.clr;
wv              = p.Results.wv;
wv_std          = p.Results.wv_std;
fs              = p.Results.fs;
saveFig         = p.Results.saveFig;

if isempty(clr)
    clr = repmat('k', 1, 2);
end

nspks = cellfun(@length, spktimes, 'uni', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gather data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ccg 50 @ 1 counts
[cc50, cc50bins] = CCG(spktimes, [], 'binSize', 0.001,...
    'duration', 0.05, 'Fs', 1 / fs);
cc50bins = cc50bins * 1000;

% ccg 100 @ 1 rate
[cc100, cc100bins] = CCG(spktimes, [], 'binSize', 0.001,...
    'duration', 0.1, 'Fs', 1 / fs);
cc100bins = cc100bins * 1000;

% ccg 20 @ 0.5 counts
[cc20, cc20bins] = CCG(spktimes, [], 'binSize', 0.0002,...
    'duration', 0.02, 'Fs', 1 / fs);
cc20bins = cc20bins * 1000;

% dccc 50 @ 1 rate
dccc = cchdeconv(cc50(:, 1, 2), cc50(:, 1, 1), nspks(1),...
    cc50(:, 2, 2), nspks(2));

% predictor
[~, pred] = cch_conv(dccc, 11, 'median', 1, 0);

% crcc
dt = diff(cc50bins(1 : 2)) / 1000;  % [s]               
den = (ones(length(cc50bins), 1) * nspks(1)) * dt; 
crccg = (dccc - pred) ./ den;

% calculate gain
roiMS = [2, 10]; 
roi_t = cc50bins >= roiMS(1) & cc50bins <= roiMS(2);
roi_t_idx = find(roi_t);
[g1, g2] = calc_stg(crccg, roi_t, dt, 0, [1, 0]);

% determine if any bin in the ROI is significant
alfa        = 0.001;
nBonf       = sum(roi_t);
gbUpper     = poissinv(1 - alfa / nBonf, max(pred(roi_t, :), [], 1));
gbLower     = poissinv(alfa / nBonf, min( pred( roi_t, : ), [], 1));
act         = any(dccc(roi_t, :) >= ones(sum(roi_t), 1) * gbUpper, 1);
sil         = any(dccc(roi_t, :) <= ones(sum(roi_t), 1) * gbLower, 1);
exc_bins    = dccc(roi_t, :) >= ones(sum(roi_t), 1) * gbUpper;
inh_bins    = dccc(roi_t, :) <= ones(sum(roi_t), 1) * gbLower;

% select excitatory or inhibitory
if sil
    stg = g2;
else
    stg = g1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMatlabGraphics(false)
fh = figure;

sb1 = subplot(2, 4, 1);     % acg100 unit1
sb2 = subplot(2, 4, 2);     % cc50 counts
sb3 = subplot(2, 4, 3);     % dccc
sb4 = subplot(2, 4, 4);     % acg100 unit2
sb5 = subplot(2, 4, 5);     % wv unit1
sb6 = subplot(2, 4, 6);     % cc100
sb7 = subplot(2, 4, 7);     % ccg25
sb8 = subplot(2, 4, 8);     % wv unit2

% acg 1
set(gcf, 'CurrentAxes', sb1)
plot_ccg(cc100(:, 1, 1), cc100bins, clr(1), [], [], [])
title(sprintf('Unit #%d (Presynaptic)', units(1)))

% acg 2
set(gcf, 'CurrentAxes', sb4)
plot_ccg(cc100(:, 2, 2), cc100bins, clr(2), [], [], [])
title(sprintf('Unit #%d (Postsynaptic)', units(2)))

% cc50 counts
set(gcf, 'CurrentAxes', sb2)
plot_ccg(cc50(:, 1, 2), cc50bins, 'k', [],...
    roi_t_idx(exc_bins), roi_t_idx(inh_bins))

% dccc
set(gcf, 'CurrentAxes', sb3)
plot_ccg(dccc, cc50bins, 'k', pred,...
    roi_t_idx(exc_bins), roi_t_idx(inh_bins))
title(sprintf('STG = %.4f', stg))

% cc100 counts 
set(gcf, 'CurrentAxes', sb6)
plot_ccg(cc100(:, 1, 2), cc100bins, 'k', [], [], [])

% cc20 counts
set(gcf, 'CurrentAxes', sb7)
plot_ccg(cc20(:, 1, 2), cc20bins, 'k', [], [], [])

if ~isempty(wv)
    
    % waveform 1
    set(gcf, 'CurrentAxes', sb5)
    x_val = [1 : size(wv, 2)] / fs * 1000;
    plot(x_val, wv(1, :), clr(1), 'LineWidth', 2)
    if ~isempty(wv_std)
        patch([x_val, flip(x_val)], [wv(1, :) + wv_std(1, :), flip(wv(1, :) - wv_std(1, :))],...
            clr(1), 'EdgeColor', 'none', 'FaceAlpha', .2, 'HitTest', 'off')
    end
    xlabel('Time [ms]')
    ylabel('Voltage [mV]')
    
    % waveform 2
    set(gcf, 'CurrentAxes', sb8)
    x_val = [1 : size(wv, 2)] / fs * 1000;
    plot(x_val, wv(2, :), clr(2), 'LineWidth', 2)
    if ~isempty(wv_std)
        patch([x_val, flip(x_val)], [wv(2, :) + wv_std(2, :), flip(wv(2, :) - wv_std(2, :))],...
            clr(2), 'EdgeColor', 'none', 'FaceAlpha', .2, 'HitTest', 'off')
    end
    xlabel('Time [ms]')
    ylabel('Voltage [mV]')
    
end
drawnow

% save
if saveFig
    figpath = fullfile(basepath, 'graphics', 'monoSyn');
    figname = fullfile(figpath, sprintf('monoSyn_%d_%d', units));
    export_fig(figname, '-jpg', '-transparent', '-r300')
end

end

% EOF

