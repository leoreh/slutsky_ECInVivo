function [bs] = getBS(varargin)

% detects burst events from LFP.
% 
% INPUT
%   sig         signal for detection
%   fs          sampling frequency
%   binsize     scalar {1}. in [s]
%   clustmet    string. cluster method. {'gmm'} or 'kmeans'
%   vars        cell of strings. one string per variable to use for
%               clustering. possible variables are 'std', 'rms', 'max',
%               'sum', 'p'
%   lognorm     logical {true}. log10 variables
%   basepath    recording session path {pwd}
%   basename    string. if empty extracted from basepath
%   graphics    logical. plot figure {1}
%   saveVar     logical. save variable {1}
%   saveFig     logical. save figure {1}
%   forceA      logical {false}. force analysis even if .mat exists
%
% OUTPUT
%   events         struct with fields:
%       stamps      n x 2 mat where n is the number of events and the
%                   columns represent the start and end of each event [samps]
%       peaks       the timestamp for each peak [samps]
%       peakPower   the voltage at each peak [uV]
%       <params>    as in input
%
% TO DO LIST
%       # GMM clusters such that grp #1 has higher values (= bursts).
%       # check rubostness binsize and interDur. Currently
%       interDur must be greater than binsize because of merge.
%       # replace changepts with https://www.mathworks.com/help/wavelet/ug/wavelet-changepoint-detection.html
%
% CALLS
%       cluDist     for cluster separation
%       calcFR
% 
% EXAMPLE
%               bs = getBS('sig', sig, 'fs', fs, 'basepath', basepath,...
%               'graphics', true, 'saveVar', false, 'binsize', 1,...
%               'clustmet', 'gmm', 'vars', {'std', 'sum', 'max'}),...
%               'saveFig', false;
%
% 31 dec 19 LH.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'sig', [], @isnumeric)
addParameter(p, 'fs', 1250, @isnumeric)
addParameter(p, 'binsize', 1, @isnumeric)
addParameter(p, 'clustmet', 'gmm', @isstr)
addParameter(p, 'vars', {'std', 'sum', 'max'}, @iscell)
addParameter(p, 'lognorm', true, @islogical)
addParameter(p, 'basepath', pwd, @isstr);
addParameter(p, 'basename', [], @isstr);
addParameter(p, 'graphics', true, @islogical)
addParameter(p, 'saveVar', true, @islogical);
addParameter(p, 'saveFig', true, @islogical);
addParameter(p, 'forceA', false, @islogical);

parse(p, varargin{:})
sig = p.Results.sig;
fs = p.Results.fs;
binsize = p.Results.binsize;
clustmet = p.Results.clustmet;
vars = p.Results.vars;
lognorm = p.Results.lognorm;
basepath = p.Results.basepath;
basename = p.Results.basename;
graphics = p.Results.graphics;
saveVar = p.Results.saveVar;
saveFig = p.Results.saveFig;
forceA = p.Results.forceAnalyze;

% params
minDur = 1;                         % minimum event duration [s]
maxDur = 300;                       % maximum event duration [s]
interDur = 1;                       % minimum time between events [bins]
binsize = binsize * fs;             % [s] * [fs] = [samples]
nclust = 2;                         % two clusters: burst and suppression

% initialize output
bs = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if mat already exists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(basename)
    [~, basename] = fileparts(basepath);
end
filename = [basepath, '\', basename, '.bs.mat'];
if ~forceA
    if exist(filename)
        load(filename)
        fprintf('\n loading %s \n\n', filename)
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% index stamps for bins
ibins = [1 : binsize : length(sig) + (mod(length(sig), binsize))];
ibins(end) = length(sig);

% divide signal to bins
sigre = sig(1 : end - (mod(length(sig), binsize) + binsize));
sigmat = reshape(sigre, binsize, (floor(length(sig) / binsize) - 1));

% last bin
siglastbin = sig(length(sigre) + 1 : length(sig));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc vars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(strcmp(vars, 'std'))
    stdvec = std(sigmat);
    stdvec = [stdvec, std(siglastbin)];
else; stdvec = []; end
if any(strcmp(vars, 'max'))
    maxvec = max(abs(sigmat));
    maxvec = [maxvec, max(abs(siglastbin))];
else; maxvec = []; end
if any(strcmp(vars, 'sum'))
    sumvec = sum(abs(sigmat));
    sumvec = [sumvec, sum(abs(siglastbin))];
else; sumvec = []; end
if any(strcmp(vars, 'rms'))
    rmsvec = rms(sigmat);
    rmsvec = [rmsvec, rms(siglastbin)];
else; rmsvec = []; end
if any(strcmp(vars, 'p'))
    [pvec, fwelch] = sum(pwelch(sigmat, 64, [], freq, 1250));
    pvec = [pvec, sum(pwelch(siglastbin, 64, [], freq, 1250))];
else; pvec = []; end

% moving
% vecwin = round(2 * fs);
% stdvec = movstd(sig, vecwin);
% sumvec = movsum(abs(sig), vecwin);
% maxvec = movmax(abs(sig), vecwin);

% concatenate and lognorm
varsmat = [stdvec; maxvec; sumvec; rmsvec; pvec];
if size(varsmat, 1) < size (varsmat, 2)
    varsmat = varsmat';
end
if lognorm
    varsmat = log10(varsmat);
end

switch clustmet
    case 'gmm'
        % fit distribution
        options = statset('MaxIter', 500);
        bs.gm = fitgmdist(varsmat, 2,...
            'options', options, 'Start', 'plus', 'Replicates', 50);
        % cluster
        [gi, ~, ~, ~, mDist] = cluster(bs.gm, varsmat);
        % cluster separation
        for i = 1 : nclust
            [lRat(i), iDist(i)] = cluDist(varsmat, find(gi == i), mDist(:, i));
        end
        
    case 'kmeans'
        [gi, cent] = kmeans(varsmat(:, varidx), 2, 'Replicates', 50, 'MaxIter', 1000);
end

% temporary fix for the arbiturary arrangment of components by gmm. assumes
% that burst has higher values for all varaibles and work for two
% components only
if mode(bs.gm.mu(1, :) > bs.gm.mu(2, :))
    gi(gi == 2) = 0;
    gi(gi == 1) = 2;
    gi(gi == 0) = 1;
    mDist = fliplr(mDist);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find start\stop times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start = find(diff(gi) > 0);
stop = find(diff(gi) < 0);

% exclude last event if it is incomplete
if length(stop) == length(start) - 1
    start = start(1 : end - 1);
end
% exclude first event if it is incomplete
if length(stop) - 1 == length(start)
    stop = stop(2 : end);
end
% correct special case when both first and last events are incomplete
if start(1) > stop(1)
    stop(1) = [];
    start(end) = [];
end

temp = [start, stop];
nevents = length(temp);

if isempty(temp)
    disp('Detection by clustering failed');
    return
else
    disp(['After clustering: ' num2str(nevents) ' events.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge events if inter-event period is too short
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iei = temp(2 : end, 1) - temp(1 : end - 1, 2);      % inter-event interval
while min(iei) < interDur
    temp2 = [];
    event = temp(1, :);
    for i = 2 : nevents
        if temp(i, 1) - event(2) < interDur
            event = [event(1) temp(i, 2)];
        else
            temp2 = [temp2; event];
            event = temp(i, :);
        end
    end
    temp2 = [temp2; event];
    temp = temp2;
    iei = temp(2 : end, 1) - temp(1 : end - 1, 2);
end

nevents = length(temp);

if isempty(temp)
    disp('Event merge failed');
    return
else
    disp(['After merging: ' num2str(nevents) ' events.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find precise sample of transition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

res = zeros(size(temp));
ipt = zeros(size(temp));
stamps = zeros(size(temp));

for i = 1 : nevents
    % start
    x = sig(ibins(temp(i, 1)) : ibins(temp(i, 1) + 2));
    [ipt(i, 1), res(i, 1)] = findchangepts(x, 'Statistic', 'std');
    temp(i, 1) = ibins(temp(i, 1)) + ipt(i, 1);
    % stop
    x = sig(ibins(temp(i, 2)) : ibins(temp(i, 2) + 2));
    [ipt(i, 2), res(i, 2)] = findchangepts(x, 'Statistic', 'std');
    temp(i, 2) = ibins(temp(i, 2)) + ipt(i, 2);
end
% from now on temp is in samples rather than bins

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discard events by duration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dur = temp(:, 2) - temp(:, 1);
idxMax = dur > maxDur * fs;
idxMin = dur < minDur * fs;
idx = idxMax | idxMin;
temp(idx, :) = [];
res(idx, :) = [];

nevents = length(temp);
if isempty(temp)
    disp('duration test failed.');
    return
else
    disp(['After testing duration: ' num2str(nevents) ' events.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% peak power and position
peakPower = zeros(nevents, 1);
peakPos = zeros(nevents, 1);
for i = 1 : nevents
    [peakPower(i), peakPos(i)] = max(abs(sig([temp(i, 1) : temp(i, 2)])));
    peakPos(i) = peakPos(i) + temp(i, 1) - 1;
end

% duration and inter-event interval
dur = temp(:, 2) - temp(:, 1);
iinterval = [temp(2 : end, 1) - temp(1 : end - 1, 2)];

% burst suppression ratio
bs.binary = zeros(length(sig), 1);
for i = 1 : size(temp, 1)
    bs.binary(temp(i, 1) : temp(i, 2)) = 1;
end
btimes = (find(~bs.binary));
bsrbinsize = 60 * fs;    % binsize for bsr calc [samples]
[bs.bsr, bs.edges, bs.cents] = calcFR(btimes, 'winCalc', [1, length(bs.binary)],...
    'binsize', bsrbinsize, 'smet', 'MA');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange struct output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bs.stamps = temp;
bs.peakPos = peakPos;
bs.peakPower = peakPower;
bs.dur = dur;
bs.iinterval = iinterval;
bs.fs = fs;
bs.binsize = binsize;
bs.vars = vars;
bs.lRat = lRat;
bs.iDist = iDist;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveVar
    save(filename, 'bs')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    f = figure;
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    suptitle(basename)
    
    % raw and BS - entire recording
    subplot(3, 4, [1 : 2])
    plot((1 : length(sig)) / fs / 60, sig)
    hold on
    axis tight
    Y = ylim;
    fill([bs.stamps fliplr(bs.stamps)] / fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
        'k', 'FaceAlpha', 0.4,  'EdgeAlpha', 0);
    ylabel('Voltage [mV]')
    set(gca, 'TickLength', [0 0])
    box off
    title('Raw and BS')
    
    % raw and BS - 1 min
    midsig = round(length(sig) / 2 / fs / 60);
    ylimsig = sig(midsig * fs * 60 : (midsig + 1) * fs * 60);
    axes('Position',[.13 .875 .15 .05])
    box on
    plot((1 : length(sig)) / fs / 60, sig)
    hold on
    axis tight
    ylim([min(ylimsig) max(ylimsig)])
    Y = ylim;
    fill([bs.stamps fliplr(bs.stamps)] / fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
        'k', 'FaceAlpha', 0.4,  'EdgeAlpha', 0);
    xlim([midsig midsig + 1])
    set(gca, 'TickLength', [0 0])
    set(gca, 'YTickLabel', [])
    set(gca, 'XTickLabel', [])
    box off
    
    % burst suppression ratio
    subplot(3, 4, [5 : 6])
    plot(bs.cents / fs / 60, bs.bsr, 'k', 'LineWidth', 3)
    ylabel('BSR [% / m]')
    xlabel('time [m]')
    axis tight
    ylim([0 1])
    set(gca, 'TickLength', [0 0])
    box off
    title('Burst Suppression Ratio')
    
    % clustering - first two variables
    options = statset('MaxIter', 50);
    gm = fitgmdist(varsmat(:, [1, 2]), 2,...
        'options', options, 'Start', 'plus', 'Replicates', 1);
    subplot(3, 4, 3)
    gscatter(varsmat(:, 1), varsmat(:, 2), gi, 'rk', '.', 1);
    axis tight
    hold on
    if strcmp(clustmet, 'gmm')
        gmPDF = @(x1, x2)reshape(pdf(gm, [x1(:) x2(:)]), size(x1));
        ax = gca;
        fcontour(gmPDF, [ax.XLim ax.YLim], 'HandleVisibility', 'off')
    end
    xlabel([vars{1} ' [log10(\sigma)]'])
    ylabel([vars{2} ' [log10(V)]'])
    legend off
    set(gca, 'TickLength', [0 0])
    box off
    title('Clustering')
    
    % clustering - second two variables
    if length(vars) > 2
        subplot(3, 4, 4)
        gm = fitgmdist(varsmat(:, [2, 3]), 2,...
            'options', options, 'Start', 'plus', 'Replicates', 1);
        gscatter(varsmat(:, 2), varsmat(:, 3), gi, 'rk', '.', 1);
        axis tight
        hold on
        if strcmp(clustmet, 'gmm')
            gmPDF = @(x1, x2)reshape(pdf(gm, [x1(:) x2(:)]), size(x1));
            ax = gca;
            fcontour(gmPDF, [ax.XLim ax.YLim], 'HandleVisibility', 'off')
        end
        xlabel([vars{2} ' [log10(V)]'])
        ylabel([vars{3} ' [log10(V)]'])
        legend({'Suppression', 'Burst'})
        set(gca, 'TickLength', [0 0])
        box off
        title(sprintf('sum(iDist) = %.2f, sum(lRat) = %.2f',...
            sum(bs.iDist), sum(bs.lRat)))
    end
    
    % power and duration with time
    subplot(3, 4, 11)
    s1 = scatter(bs.peakPos / fs / 60, bs.dur / fs, 2, 'k', 'filled');
    l1 = lsline;
    ylabel('Duration [s]')
    yyaxis right
    s2 = scatter(bs.peakPos / fs / 60, bs.peakPower, 2, 'b', 'filled');
    l2 = lsline;
    set(l1, 'color', 'k')
    set(l2, 'color', 'b')
    l1.LineWidth = 3;
    l2.LineWidth = 3;
    axis tight
    ylabel('Peak Voltage [uV]')
    xlabel('Time [m]')
    set(gca, 'TickLength', [0 0])
    box off
    title('Parms vs. Time')

    % duration and power vs inter-burst interval
    subplot(3, 4, 12)
    s1 = scatter((bs.iinterval) / fs, bs.dur(2 : end) / fs, 2, 'k', 'filled');
    l1 = lsline;
    yyaxis right
    s2 = scatter(bs.iinterval / fs, bs.peakPower(2 : end), 2, 'b', 'filled');
    l2 = lsline;
    set(l1, 'color', 'k')
    set(l2, 'color', 'b')
    l1.LineWidth = 3;
    l2.LineWidth = 3;
    set(gca, 'TickLength', [0 0])
    box off
    axis tight
    xlim([0 30])
    ylabel('Duration [s]')
    xlabel('Interval [s]')
    title('Inter-burst interval')
    
    % duration histogram
    subplot(3, 4, 7)
    h = histogram(log10(bs.dur / fs), 30, 'Normalization', 'Probability');
    h.EdgeColor = 'none';
    h.FaceColor = 'k';
    h.FaceAlpha = 0.5;
    hold on
    plot(log10([minDur minDur]), ylim, '--k')
    xlabel('Time [log(s)]')
    ylabel('Probability [%]')
    set(gca, 'TickLength', [0 0])
    box off
    title('Burst duration')
    
    % amplitude histogram
    subplot(3, 4, 8)
    h = histogram(log10(bs.peakPower), 50, 'Normalization', 'Probability');
    h.EdgeColor = 'none';
    h.FaceColor = 'k';
    h.FaceAlpha = 0.5;
    xlabel('Peak voltage [log(uV)]')
    ylabel('Probability [%]')
    set(gca, 'TickLength', [0 0])
    box off
    title('Burst amplitude')
    
    if saveFig
        figname = [basename '_BS'];
        export_fig(figname, '-tif', '-transparent')
        % savePdf(figname, basepath, ff)
    end
end
end

% EOF

