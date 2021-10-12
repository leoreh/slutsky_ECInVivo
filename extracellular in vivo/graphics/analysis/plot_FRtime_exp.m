
% plots firing rate of mu (sr) and su (pyr and int) across time
% concatenated from multiple sessions. depends on varArray.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grp = [1 : 4];                  % which tetrodes to plot
suFlag = 1;                     % plot only su or all units
frBoundries = [0 Inf];          % include only units with fr greater / lower than
saveFig = false;
sessionIdx = 1 : nsessions;     % selection of sessions

% -------------------------------------------------------------------------

% find date time of start and end of experiment (round to hour)
assignVars(varArray, sessionIdx(end))
ts = sr.binsize;     
[dtStart, ~] = guessDateTime(dirnames(sessionIdx(1)));
dtStart = dateshift(dtStart, 'start', 'hour');
dtStart = dtStart - seconds(60 * 60); % override
[dtEnd, ~] = guessDateTime(dirnames(sessionIdx(end)));
dtEnd = dtEnd + seconds(length(sr.strd) * ts);
dtEnd = dateshift(dtEnd, 'end', 'hour');
expLen = seconds(dtEnd - dtStart) / ts;     

% get datetime indices and labels for x axis
dtAxis = dtStart : hours(6) : dtEnd;
tidx = 1;
for itime = 1 : length(dtAxis)
    [~, tidx(itime)] = tstamp2time('dtstr', dtStart, 'tstr', dtAxis(itime), 'fs', 1 / ts);
    if tidx(itime) == 0
        tidx(itime) = 1;
    end
    tlabel{itime} = datestr(datenum(dtAxis(itime)), 'dd/mm_HH:MM');
end
 
% get tidx of manipulation
dt1 = datetime(2021, 02, 28, 19, 00, 00);
dt2 = datetime(2021, 03, 03, 17, 10, 00);
[~, shadeIdx(1)] = tstamp2time('dtstr', dtStart, 'tstr', dt1, 'fs', 1 / ts);
[~, shadeIdx(2)] = tstamp2time('dtstr', dtStart, 'tstr', dt2, 'fs', 1 / ts);
shadeIdx = [0, 0];  % override

% initialize vars that will carry information from entire experiment
clear units
for isession = 1 : nsessions
    assignVars(varArray, isession)
    units(isession, 1) = sum(selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'pyr'));
    units(isession, 2) = sum(selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'int'));
end
maxUnitsSession = max(units);
expRS = nan(expLen, maxUnitsSession(1));
expFS = nan(expLen, maxUnitsSession(2));
expMU = nan(expLen, length(grp));

% concatenate data from different sessions. finds the index of rec start
% from basename and assumes the session recording is continuous
idx_recStart = zeros(length(sessionIdx), 1);
for isession = sessionIdx
    assignVars(varArray, isession)   
    
    % initialize
    RSunits{isession} = nan(1, 1);
    FSunits{isession} = nan(1, 1);

    % indices
    idx_recStart(isession) = floor(max([1, seconds(guessDateTime(dirnames(isession)) - dtStart) / ts]));
    idx_session = idx_recStart(isession) : idx_recStart(isession) + length(sr.strd) - 1;
    
    % MU
    expMU(idx_session, :) = sr.strd(grp, :)';

    if ~isempty(varArray{isession, 3})
        % RS
        RSunits{isession} = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'pyr');
        expRS(idx_session, 1 : sum(RSunits{isession})) = fr.strd(RSunits{isession}, :)';
        
        % FS
        FSunits{isession} = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'int');
        expFS(idx_session, 1 : sum(FSunits{isession})) = fr.strd(FSunits{isession}, :)';
    end 
end

% plot
fh = figure;
sb1 = subplot(2, 1, 1);
yLimit = [0 max(max(expMU))];
rectangle('Position', [shadeIdx(1) yLimit(1) diff(shadeIdx) diff(yLimit)],...
    'Curvature',0.2, 'faceColor', [0.95 0.95 0.95])
hold on
plot(expMU)
axis tight
xlim([1 length(expMU)])
set(gca, 'box', 'off');
xticks(tidx)
xticklabels(tlabel)
xtickangle(45)
ylabel('Multi-unit firing rate [Hz]')
lgh = legend(split(num2str(grp)));
ylim([0 100])

sb2 = subplot(2, 1, 2);
plot(mean(expRS', 'omitnan'), 'b', 'LineWidth', 2)
hold on
plot(mean(expFS', 'omitnan'), 'r', 'LineWidth', 2)
yLimit = [0 max(mean(expFS', 'omitnan'))];
plot([idx_recStart, idx_recStart], yLimit, '--k');
axis tight
xlabel('Time [h]')
ylabel('Single unit firing rate [Hz]')
set(gca, 'box', 'off')
xticks(tidx)
xticklabels(tlabel)
xtickangle(45)
linkaxes([sb1, sb2], 'x')
legend(sprintf('RS ~= %d su', round(mean(sum(cell2nanmat(RSunits(sessionIdx), 2), 'omitnan')))),...
 sprintf('FS ~= %d su', round(mean(sum(cell2nanmat(FSunits(sessionIdx), 2), 'omitnan')))));
if saveFig
    figpath = fullfile(mousepath, 'graphics');
    mkdir(figpath)
    figname = fullfile(figpath, ['FR_time']);
    export_fig(figname, '-tif', '-transparent', '-r300')
end

% -------------------------------------------------------------------------
%%% find tidx from sample in specific dat file
% sessionIdx2 = 1;
% basename = dirnames(sessionIdx2);
% assignVars(varArray, sessionIdx2)
% csum = cumsum(datInfo.nsamps);
% [dtSession, ~] = guessDateTime(basename);
% 
% % get time within session
% [dt, ~] = tstamp2time('dtstr', basename, 'tstamp', csum(2), 'fs', fs);
% % get idx relative to experiment start
% idx = (hours(dt - dtStart))';
% %%%
% 
% %%% x axis in hr
% xAxisHr = [1 : length(expFS)] / ts;
% %%%