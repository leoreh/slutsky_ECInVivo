function plot_monoSyn(varargin)

% plot mono synaptic connection
% 
% INPUT:
%   basepath        path to recording {pwd}
%   spktimes        cell array of spike times [s]
%   units           1 x 2 numeric. indices to spktimes
%   clr             2 element char. color for each cell
%   swv             struct. see spkwvMetrics.m
%   fs              numeric. sampling frequency 
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
addOptional(p, 'units', [1, 2], @isnumeric);
addOptional(p, 'clr', 'kk', @ischar);
addOptional(p, 'swv', []);
addOptional(p, 'fs', 10000, @isnumeric);

parse(p, varargin{:})
basepath        = p.Results.basepath;
spktimes        = p.Results.spktimes;
units           = p.Results.units;
clr             = p.Results.clr;
swv             = p.Results.swv;
fs              = p.Results.fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gather data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(clr)
    clr = repmat('k', 1, 2);
end

acg_wide_bnsz = 0.0001;
acg_wide_dur = 0.02;
acg_narrow_bnsz = 0.0005;
acg_narrow_dur = 0.1;

% acg narrow
[acg_narrow, acg_narrow_tstamps] = CCG(spktimes(units),...
    [], 'binSize', acg_narrow_bnsz,...
    'duration', acg_narrow_dur, 'norm', 'rate', 'Fs', 1 / fs);

% ccg super narrow
[acg_wide, acg_wide_tstamps] = CCG(spktimes(units),...
    [], 'binSize', acg_wide_bnsz,...
    'duration', acg_wide_dur, 'norm', 'rate', 'Fs', 1 / fs);

acg_wide_tstamps = acg_wide_tstamps * 1000;
acg_narrow_tstamps = acg_narrow_tstamps * 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMatlabGraphics(false)
fh = figure;

% acg narrow 1
subplot(2, 3, 1)
bh = bar(acg_narrow_tstamps,...
    squeeze(acg_narrow(:, 1, 1)), 'BarWidth', 1);
ylabel('Rate [Hz]')
xlabel('Time [ms]')
bh.FaceColor = clr(1);
bh.EdgeColor = 'none';
box off
axis tight
subtitle(sprintf('Unit #%d (Presynaptic)', units(1)))

% ccg 
subplot(2, 3, 2)
bh = bar(acg_narrow_tstamps,...
    squeeze(acg_narrow(:, 1, 2)), 'BarWidth', 1);
hold on
plot([0, 0], ylim, '--k')
ylabel('Rate [Hz]')
xlabel('Time [ms]')
bh.FaceColor = 'k';
bh.EdgeColor = 'none';
box off
axis tight

% acg narrow 2
subplot(2, 3, 3)
bh = bar(acg_narrow_tstamps,...
    squeeze(acg_narrow(:, 2, 2)), 'BarWidth', 1);
ylabel('Rate [Hz]')
xlabel('Time [ms]')
bh.FaceColor = clr(2);
bh.EdgeColor = 'none';
box off
axis tight
subtitle(sprintf('Unit #%d (Postsynaptic)', units(2)))

% ccg super narrow 
subplot(2, 3, 5)
bh = bar(acg_wide_tstamps,...
    squeeze(acg_wide(:, 1, 2)), 'BarWidth', 1);
hold on
plot([0, 0], ylim, '--k')
ylabel('Rate [Hz]')
xlabel('Time [ms]')
bh.FaceColor = 'k';
bh.EdgeColor = 'none';
box off
axis tight

if ~isempty(swv)
    
    % waveform 1
    subplot(2, 3, 4)
    x_val = [1 : size(swv.wv, 2)] / fs * 1000;
    wv = swv.wv(units(1), :);
    wv_std = swv.wv_std(units(1), :);
    plot(x_val, wv, clr(1), 'LineWidth', 2)
    patch([x_val, flip(x_val)], [wv + wv_std, flip(wv - wv_std)],...
        clr(1), 'EdgeColor', 'none', 'FaceAlpha', .2, 'HitTest', 'off')
    xlabel('Time [ms]')
    ylabel('Voltage [mV]')
    
    % waveform 2
    subplot(2, 3, 6)
    x_val = [1 : size(swv.wv, 2)] / fs * 1000;
    wv = swv.wv(units(2), :);
    wv_std = swv.wv_std(units(2), :);
    plot(x_val, wv, clr(2), 'LineWidth', 2)
    patch([x_val, flip(x_val)], [wv + wv_std, flip(wv - wv_std)],...
        clr(2), 'EdgeColor', 'none', 'FaceAlpha', .2, 'HitTest', 'off')
    xlabel('Time [ms]')
    ylabel('Voltage [mV]')
    
    % save
    figpath = fullfile(basepath, 'graphics');
    figname = fullfile(figpath, sprintf('monoSyn_%d_%d', units));
    export_fig(figname, '-jpg', '-transparent', '-r300')
    
end

end

% EOF