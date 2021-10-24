
% plots firing rate of mu (sr) and su (pyr and int) across time for a
% single session (according to basepath).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = pwd;
[~, basename] = fileparts(basepath);
cd(basepath)

grp = [1 : 4];                  % which tetrodes to plot
suFlag = 1;                     % plot only su or all units
frBoundries = [0 Inf];          % include only units with fr greater / lower than
saveFig = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vars = ["session.mat";...
    "cell_metrics.cellinfo";...
    "spikes.cellinfo";...
    "fr.mat";...
    "datInfo";...
    "AccuSleep_states";...
    "sr.mat"];

[varArray, ~, mousepath] = getSessionVars('vars', vars,...
    'dirnames', string(basename));
assignVars(varArray, 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(fr.strd) ~= length(sr.strd)
    warning('check length of mu and su firing rate')
end

% x axis in hr
ts = fr.binsize;
xidx = [1 : length(fr.strd)] / ts;

% idx of block tranisition (dashed lines)
if ~isempty(datInfo)
    if ~isfield('datInfo', 'fs')
        datInfo.fs = 20000;
    end
    csum = cumsum(datInfo.nsamps) / datInfo.fs / 60 / 60;
    tidx = csum(1);
else
    tidx = 0;
end

% units
clear units
units(1, :) = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'pyr');
units(2, :) = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'int');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure;
sb1 = subplot(2, 1, 1);
yLimit = ceil([0 max(max(sr.strd))]);
hold on
plot(xidx, sr.strd')
plot([tidx tidx], yLimit, '--k')
axis tight
set(gca, 'box', 'off');
ylabel('Multi-unit firing rate [Hz]')
lgh = legend(split(num2str(grp)));

sb2 = subplot(2, 1, 2);
yLimit = [0 ceil(max(mean(fr.strd(units(2, :), :), 'omitnan')))];
hold on
plot(xidx, mean(fr.strd(units(1, :), :), 'omitnan'), 'b', 'LineWidth', 2)
plot(xidx, mean(fr.strd(units(2, :), :), 'omitnan'), 'r', 'LineWidth', 2)
plot([tidx tidx], yLimit, '--k')
axis tight
xlabel('Time [h]')
ylabel('Single unit firing rate [Hz]')
set(gca, 'box', 'off')
linkaxes([sb1, sb2], 'x')
legend(sprintf('RS ~= %d su', sum(units(1, :))),...
    sprintf('FS ~= %d su', sum(units(2, :))));

if saveFig
    figpath = fullfile(pwd, 'graphics');
    mkdir(figpath)
    figname = fullfile(figpath, ['FR_time']);
    export_fig(figname, '-tif', '-transparent', '-r300')
end

guessDateTime(basename)
