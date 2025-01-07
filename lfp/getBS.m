function [bs] = getBS(varargin)

% receives an lfp signal, divides it to 500 ms bins and for each bin
% calculates several variables (e.g. the first principle component of the
% spectrogram). clusters each bin as either a burst or suppression. defines
% bouts of burst while accounting for duration parameters. finds the
% precise timing of a transition between an bout of burst and suppression
% (currently not implemented). allows the user to go over the data
% and manually mark bouts of burst. calculates several params such as the
% burst suppression ratio. 
% 
% INPUT
%   sig         numeric vec of lfp data.
%   fs          sampling frequency {1250}.
%   smf         smooth factor for BSR and delta power [bins] {7}.
%   binsize     scalar [s] {0.5}.
%   BSRbinsize  scalar [samples]. default {65536} to fit spectrogram 
%               w/ 10 s window and 0% overlap
%   clustmet    string. cluster method. {'gmm'} or 'kmeans'
%   vars        cell of strings. one string per variable to use for
%               clustering. possible variables are 'std', 'rms', 'max',
%               'sum', 'p', 'energy', 'pc1'
%   lognorm     logical {true}. log10 variables.
%   basepath    recording session path {pwd}
%   basename    string. if empty extracted from basepath
%   graphics    logical. plot figure {1}
%   saveVar     logical. save variable {1}
%   saveFig     logical. save figure {1}
%   vis         logical. figure visible {1} or not
%   forceA      logical {false}. force analysis even if .mat exists
%   mancur      logical {false}. manually mark bursts.
%
% OUTPUT
%   bs              struct with fields:
%       stamps      n x 2 mat where n is the number of events and the
%                   columns represent the start and end of each event [samps]
%       manStamps   same as stamps but of manually selected bursts
%       accuracy    percent accuracy of manual detection vs. algorithm
%       peaksPos    the timestamp for each peak [samps]
%       peakPower   the voltage at each peak [uV]
%       edges       of bins for bsr
%       cents       of bins for bsr
%       lRat        L-ratio of clusters separation
%       iDist       isolation distance between clusters
%       bsr         burst suppression ratio
%       dur         burst duration [samples]
%       iinterval   inter-burst interval [samples]
%       gm          gaussian mixture model used for separation
%       binary      binary vector where 1 = burst and 0 = suppression
%       <params>    as in input
%
% TO DO LIST
%       # GMM clusters such that grp #1 has higher values (= bursts).
%       # account for robustness to minDur / interDur. 
%         Currently the smallest minDur is binsize
%       # replace changepts with https://www.mathworks.com/help/wavelet/ug/wavelet-changepoint-detection.html
%
% CALLS
%       cluDist             cluster separation
%       times2rate          BSR 
%       specBand            delta band
%       binary2bouts       start/stop times
% 
% EXAMPLE
%        bs = getBS('sig', sig, 'fs', fs, 'basepath', basepath,...
%             'graphics', true, 'saveVar', true, 'binsize', 0.5,...
%             'clustmet', 'gmm', 'vars', vars, 'basename', basename,...
%             'saveFig', saveFig, 'forceA', forceA);
%
% 31 dec 19 LH  updates
% 06 may 20 LH      manual selection of bursts
% 03 aug 20 LH      times2rate instead of calcFR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'sig', [], @isnumeric)
addParameter(p, 'fs', 1250, @isnumeric)
addParameter(p, 'smf', 7, @isnumeric)
addParameter(p, 'binsize', 0.5, @isnumeric)
addParameter(p, 'BSRbinsize', 65536, @isnumeric)
addParameter(p, 'clustmet', 'gmm', @isstr)
addParameter(p, 'vars', {'std', 'pc1', 'max'}, @iscell)
addParameter(p, 'lognorm', true, @islogical)
addParameter(p, 'basepath', pwd, @isstr);
addParameter(p, 'basename', [], @isstr);
addParameter(p, 'graphics', true, @islogical)
addParameter(p, 'saveVar', true, @islogical);
addParameter(p, 'saveFig', true, @islogical);
addParameter(p, 'vis', true, @islogical);
addParameter(p, 'forceA', false, @islogical);
addParameter(p, 'mancur', false, @islogical);

parse(p, varargin{:})
sig = p.Results.sig;
fs = p.Results.fs;
smf = p.Results.smf;
binsize = p.Results.binsize;
BSRbinsize = p.Results.BSRbinsize;
clustmet = p.Results.clustmet;
vars = p.Results.vars;
lognorm = p.Results.lognorm;
basepath = p.Results.basepath;
basename = p.Results.basename;
graphics = p.Results.graphics;
saveVar = p.Results.saveVar;
saveFig = p.Results.saveFig;
vis = p.Results.vis;
forceA = p.Results.forceA;
mancur = p.Results.mancur;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minDur = 1;                             % minimum event duration [bins]
maxDur = 3000000;                       % maximum event duration [bins]
interDur = 4;                           % minimum time between events [bins]
binsize = 2 ^ nextpow2(binsize * fs);   % for pc1
nclust = 2;                             % 2 clusters: burst and suppression

% initialize output
bs.bsr = []; bs.edges = []; bs.cents = []; bs.cents = []; bs.stamps = [];
bs.manStamps = []; bs.accuracy = []; bs.peakPos = []; bs.peakPower = []; 
bs.dur = []; bs.iinterval = []; bs.fs = fs; bs.binsize = binsize;
bs.lRat = []; bs.iDist = []; bs.vars = vars;

if isempty(basename)
    [~, basename] = fileparts(basepath);
end
filename = [basepath, '\', basename, '.bs.mat'];
if ~forceA
    if exist(filename)
        load(filename)
        fprintf('\n loading %s \n', filename)
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% index stamps for bins
ibins = [1 : binsize : length(sig)];
ibins(end) = length(sig);

% divide signal to bins
sigre = sig(1 : end - (mod(length(sig), binsize) + binsize));
sigmat = reshape(sigre, binsize, (floor(length(sig) / binsize) - 1));
siglastbin = sig(length(sigre) + 1 : length(sig)); % last bin

% calc vars
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
    [pvec] = sum(abs(pwelch(sigmat, 64, [], logspace(0, 2, 100), 1250)));
    pvec = [pvec, sum(abs(pwelch(siglastbin, 64, [], logspace(0, 2, 100), 1250)))];
else; pvec = []; end
if any(strcmp(vars, 'energy'))
    for i = 1 : length(ibins) - 1
        evec(i) = sum(abs(energyop(sig(ibins(i) : ibins(i + 1)), false)));
    end
else; evec = []; end
if any(strcmp(vars, 'pc1'))
    freq = logspace(0, 2, 100);
    win = hann(binsize);
    [~, fff, ttt, pband] = spectrogram(sig, win, 0, freq, fs, 'yaxis', 'psd');
    pband = 10 * log10(abs(pband));
    [~, pc1] = pca(pband', 'NumComponents', 1);
else; evec = []; end

% concatenate and lognorm
varmat = [stdvec; maxvec; sumvec; rmsvec; pvec; evec];
if size(varmat, 1) < size (varmat, 2)
    varmat = varmat';
end
if lognorm
    varmat = log10(varmat);
end
varmat = [varmat, pc1];

switch clustmet
    case 'gmm'
        
        % fit distribution
        options = statset('MaxIter', 500);
        bs.gm = fitgmdist(varmat, nclust,...
            'options', options, 'Start', 'plus', 'Replicates', 50);
        
        % cluster
        [gi, ~, ~, ~, mDist] = cluster(bs.gm, varmat);
        
        % cluster separation
        for i = 1 : nclust
            [lRat(i), iDist(i)] = cluDist(varmat, find(gi == i), mDist(:, i));
        end
        
    case 'kmeans'
        [gi, cent] = kmeans(varmat(:, varidx), 2, 'Replicates', 50, 'MaxIter', 1000);
end

% temporary fix for the arbiturary arrangment of components by gmm. assumes
% that bursts have higher values for all varaibles. works for two
% components only. Thus grp 1 = suppression and grp 2 = brust
if mode(bs.gm.mu(1, :) > bs.gm.mu(2, :))
    gi(gi == 2) = 0;
    gi(gi == 1) = 2;
    gi(gi == 0) = 1;
    mDist = fliplr(mDist);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find start\stop times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stamps = binary2bouts('vec', gi, 'minDur', minDur, 'maxDur', maxDur,...
    'interDur', interDur);
stamps = ibins(stamps);
nevents = size(stamps, 1);

if nevents == 0
    fprintf('\n\nDetection by clustering failed\n\n');
    bs.bsr = zeros(1, floor(length(sig) / binsize));
    bs.cents = zeros(1, floor(length(sig) / binsize));
    if saveVar
        save(filename, 'bs')
    end
    return
else
    fprintf('\n\nAfter clustering: %d events\n\n', nevents);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find precise sample of transition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% res = zeros(size(stamps, 1), 1);
% for i = 1 : nevents
%     % OPTION 1: find one change point for start and stop separately
%     marg = binsize * 2;
%     % start
%     idx = max([stamps(i, 1) - marg, 1]);
%     x = sig(idx : idx + 2 * marg);
%     [ipt, res(i, 1)] = findchangepts(x, 'Statistic', 'rms');
%     if ~isempty(ipt)
%         stamps(i, 1) = idx + ipt;
%         % correct intances where findchangepts detected decreases in std
%         if rms(x(1 : ipt)) > rms(x(ipt + 1 : end))
%             stamps(i, 1) = 0;
%         end
%     end
%     % stop
%     idx = min([stamps(i, 2) + marg, length(sig)]);
%     x = sig(stamps(i, 2) - marg : idx);
%     [ipt, res(i, 2)] = findchangepts(x, 'Statistic', 'rms');
%     if ~isempty(ipt)
%         stamps(i, 2) = stamps(i, 2) - marg + ipt;
%         % correct intances where findchangepts detected increases in std
%         if rms(x(i : ipt)) < rms(x(ipt + 1 : end))
%             stamps(i, 2) = 0;
%         end
%     end
%     
% % %     % OPTION 2: find two change points per burst. very slow
% % %     marg = binsize * 2;
% % %     idx(1) = max([stamps(i, 1) - marg, 1]);
% % %     idx(2) = min([stamps(i, 2) + marg, length(sig)]);
% % %     x = sig(idx(1) : idx(2));
% % % 
% % %     [ipt, res(i)] = findchangepts(x, 'Statistic', 'std',...
% % %         'MaxNumChanges', 2, 'MinDistance', minDur * fs);  
% % % 
% % %     % remove bursts with < 2 change points.
% % %     if length(ipt) == 2
% % %         stamps(i, 1) = ibins(temp(i, 1) - 1) + ipt(1);
% % %         stamps(i, 2) = ibins(temp(i, 1) - 1) + ipt(2);
% % %         % correct intances where findchangepts detected decreases in std
% % %         if std(sig(stamps(i, 1) : stamps(i, 2))) < std(x)
% % %             stamps(i, :) = [0 0];
% % %         end
% % %     end   
% 
% end
% stamps(any(stamps' == 0), :) = [];
% res(any(stamps' == 0), :) = [];
% 
% nevents = size(stamps, 1);
% if isempty(stamps)
%     fprintf('Localization failed');
%     return
% else
%     fprintf('After localization: %d events\n', nevents);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% peak power and position
peakPower = zeros(nevents, 1);
peakPos = zeros(nevents, 1);
for i = 1 : nevents
    [peakPower(i), peakPos(i)] = max(abs(sig(stamps(i, 1) : stamps(i, 2))));
    peakPos(i) = peakPos(i) + stamps(i, 1) - 1;
end

% duration and inter-event interval
dur = stamps(:, 2) - stamps(:, 1);
iinterval = [stamps(2 : end, 1) - stamps(1 : end - 1, 2)];

% burst suppression ratio
bs.binary = zeros(length(sig), 1);
for i = 1 : size(stamps, 1)
    bs.binary(stamps(i, 1) : stamps(i, 2)) = 1;
end
btimes = (find(~bs.binary));
[bs.bsr, bs.edges, bs.cents] = times2rate(btimes, 'winCalc',...
    [1, length(bs.binary)], 'binsize', BSRbinsize, 'c2r', true);
bs.bsr = movmean(bs.bsr, smf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manual detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mancur
    
    % mark episodes of bursts
    bs.manStamps = markEp((1 : length(sig)) / fs / 60, sig);
    
    % convert episodes of bursts to binary vector
    x = floor(bs.manStamps * fs) * 60;
    epbin = zeros(length(sig), 1);
    for i = 1 : size(x, 1)
        epbin(x(i, 1) : x(i, 2)) = 1;
    end
    
    % correct special case where user marked beyond end of recording. very
    % laze fix
    epbin(length(bs.binary) + 1 : length(epbin)) = [];
    
    % sensitivity (TPR) 
    tpr = sum(~bs.binary & ~epbin) / sum(~epbin);
    
    % specificity (TNR)
    tnr = sum(bs.binary & epbin) / sum(epbin);
    
    % accuracy
    bs.accuracy = (sum(bs.binary & epbin) + sum(~bs.binary & ~epbin)) / length(sig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange struct output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bs.stamps = stamps;
bs.peakPos = peakPos;
bs.peakPower = peakPower;
bs.dur = dur;
bs.iinterval = iinterval;
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
    if vis
        fh = figure;
    else
        fh = figure('Visible', 'off');
    end
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    suptitle(basename)

    % raw and BS - entire recording
    sb1 = subplot(3, 4, [1 : 2]);
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
    
    % burst suppression ratio, delta
    sb2 = subplot(3, 4, [5 : 6]);
    plot(bs.cents / fs / 60, bs.bsr, 'k', 'LineWidth', 2)
    hold on
    % delta band power
    [dband, ~] = specBand('sig', sig, 'graphics', false,...
        'band', [1 4], 'binsize', BSRbinsize, 'smf', smf);
    plot(bs.cents / fs / 60, dband, 'r', 'LineWidth', 2)
    legend({'BSR', '[1-4 Hz]'})
    ylabel('[a.u]')
    xlabel('Time [m]')
    axis tight
    ylim([0 1])
    set(gca, 'TickLength', [0 0])
    box off
    title('Delta power and BSR')
    
    % spectrogram
    %     sb3 = subplot(3, 4, [9 : 10]);
    %     specBand('sig', sig, 'band', [], 'graphics', true,...
    %         'binsize', BSRbinsize, 'smf', smf);
    
    % clustering - first two variables
    options = statset('MaxIter', 50);
    gm = fitgmdist(varmat(:, [1, 2]), nclust,...
        'options', options, 'Start', 'plus', 'Replicates', 50);
    subplot(3, 4, 3)
    gscatter(varmat(:, 1), varmat(:, 2), gi, 'rk', '.', 1);
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
        gm = fitgmdist(varmat(:, [2, 3]), 2,...
            'options', options, 'Start', 'plus', 'Replicates', 1);
        gscatter(varmat(:, 2), varmat(:, 3), gi, 'rk', '.', 1);
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
    ylabel('Duration [s]')
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
    ylabel('Peak Voltage [uV]')
    xlabel('Interval [s]')
    title('Inter-burst interval')
    
    % duration histogram
    subplot(3, 4, 7)
    h = histogram(log10(bs.dur / fs), 30, 'Normalization', 'Probability');
    h.EdgeColor = 'none';
    h.FaceColor = 'k';
    h.FaceAlpha = 1;
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
    h.FaceAlpha = 1;
    xlabel('Peak voltage [log(uV)]')
    ylabel('Probability [%]')
    set(gca, 'TickLength', [0 0])
    box off
    title('Burst amplitude')
    
%     linkaxes([sb1, sb2, sb3], 'x');
    
    if saveFig
        figname = [basename '_BS'];
        export_fig(figname, '-tif', '-transparent')
    end
end
end

% EOF

