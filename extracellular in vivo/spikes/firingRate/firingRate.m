function fr = firingRate(spktimes, varargin)

% wrapper for firing rate functions. calculates firing rate based on spike
% times and smoothes result by a moving average (MA) or Gaussian kernel
% (GK) impleneted by multiple-pass MA. Default is to calculate firing rate
% in sliding 1-min windows of 20 s steps (Miyawaki et al., Sci. Rep.,
% 2019). In practice this is done by setting binsize to 60 s and smoothing
% w/ moving average of 3 points.
% 
% INPUT
%   spktimes    a cell array of vectors. each vector (unit) contains the
%               timestamps of spikes. for example {spikes.times{1:4}}
%   basepath    recording session path {pwd}
%   graphics    plot figure {1}.
%   saveFig     save figure {1}.
%   saveVar     save variable {1}.
%   winCalc     time window for calculation {[1 Inf]}. specified in s.
%   binsize     size bins {60}. specified in s.
%   metBL       calculate baseline as 'max' or {'avg'}.
%   winBL       window to calculate baseline FR {[1 Inf]}.
%               specified in s.
%   select      cell array with strings expressing method to select units.
%               'thr' - units with fr > 0.05 during baseline
%               'stable' - units with std of fr < avg of fr during
%               baseline. default = none.
%   smet        method for smoothing firing rate: moving average (MA) or
%               Gaussian kernel (GK) impleneted by multiple-pass MA. {[]}.
% 
% OUTPUT
% fr            struct with fields strd, norm, bins, binsize,
%               normMethod, normWin
%
% 26 feb 19 LH. 
% 21 nov 19 LH  added active periods and mFR accordingly
% 09 may 20 LH  fixed c2r issues with calcFR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validate_win = @(win) assert(isnumeric(win) && length(win) == 2,...
    'time window must be in the format [start end]');

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'binsize', 60, @isscalar);
addOptional(p, 'winCalc', [1 Inf], validate_win);
addOptional(p, 'winBL', [], validate_win);
addOptional(p, 'metBL', 'avg', @ischar);
addOptional(p, 'select', {''}, @iscell);
addOptional(p, 'smet', 'none', @ischar);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveFig', true, @islogical);
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
binsize = p.Results.binsize;
winCalc = p.Results.winCalc;
winBL = p.Results.winBL;
metBL = p.Results.metBL;
select = p.Results.select;
smet = p.Results.smet;
graphics = p.Results.graphics;
saveFig = p.Results.saveFig;
saveVar = p.Results.saveVar;

% validate window
if winCalc(2) == Inf
    winCalc(2) = max(vertcat(spktimes{:}));
end

smfactor = 9;    % smooth factor
nunits = length(spktimes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate firing rate
[fr.strd, ~, fr.tstamps] = times2rate(spktimes, 'binsize', binsize,...
    'winCalc', winCalc, 'c2r', true);

% smooth
switch smet
    case 'MA'
        fr.strd = movmean(fr.strd, smfactor, 2);
    case 'GK'
        gk = gausswin(smfactor);
        gk = gk / sum(gk);
        for i = 1 : nunits
            x(i, :) = conv(fr.strd(i, :), gk, 'same');
        end
end

% normalize firing rate
if isempty(winBL)
    winBL = [1 size(fr.strd, 2)];
else
    winBL = winBL / binsize;
    if winBL(1) < 1
        winBL(1) = 1;
    end
    if winBL(2) == Inf
        winBL(2) = size(fr.strd, 2);
    end
end
fr.norm = fr.strd ./ mean(fr.strd(:, winBL(1) : winBL(2)), 2);

% select units who fired above thr
if any(strcmp(select, 'thr'))
    ithr = bl > 0.05;
else
    ithr = ones(nunits, 1);
end
% select units with low variability
if any(strcmp(select, 'stable'))
    bl_std = std(fr(:, win(1) : win(2)), [], 2);
    istable = bl_std < bl;
else
    istable = ones(nunits, 1);
end
idx = istable & ithr;

% add params
fr.winBL = winBL;
fr.binsize = binsize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphics
    plotFRtime('fr', fr, 'spktimes', spktimes, 'units', true,...
        'avg', true, 'raster', true, 'saveFig', saveFig)  
    
%     bl = blFR(fr.norm, 'method', metBL, 'win', winBL);
%     plotFRdistribution(bl, 'saveFig', saveFig) 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveVar  
    [~, filename] = fileparts(basepath);
    save([basepath, filesep, filename, '.fr.mat'], 'fr')
end

return

% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% messaround with figuress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fs = 20000;
fs = 24414.1;
binsize = 60;

% select units
pyr = strcmp(cell_metrics.putativeCellType, 'Pyramidal Cell');
int = strcmp(cell_metrics.putativeCellType, 'Narrow Interneuron');
su = spikes.su';
grp = [1 : 4];
grpidx = zeros(1, length(spikes.shankID));
for ii = 1 : length(grp)
    grpidx = grpidx | spikes.shankID == grp(ii);
end

% override 
su = ones(1, length(spikes.ts));
pyr = ones(1, length(spikes.su));

% states
states = {SleepState.ints.WAKEstate, SleepState.ints.NREMstate, SleepState.ints.REMstate};
for ii = 1 : length(states)
    tStates{ii} = InIntervals(fr.tstamps, states{ii});
    t{ii} = fr.tstamps(tStates{ii});
    frStates{ii} = mean(fr.strd(:, tStates{ii}), 2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare two time periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
units = pyr & su & grpidx;

win = floor([1, datInfo.nsamps(1) / fs; datInfo.nsamps(1) / fs sum(datInfo.nsamps) / fs]);
win(1) = 1;
win(2) = win(2) + 1;

wint1 = InIntervals(fr.tstamps, win(1, :));
wint2 = InIntervals(fr.tstamps, win(2, :));

idx1 = wint1 & tStates{1};
idx2 = wint2 & tStates{1};

avg1 = mean(fr.strd(units, idx1), 2);
avg2 = mean(fr.strd(units, idx2), 2);

% all units
avg1 = mean(fr.strd(:, :), 2);
avg2 = mean(fr.strd(:, :), 2);


c = repmat([0 0 1], length(spikes.su), 1);
c(int, :) = repmat([1 0 0], sum(int), 1);
scatter(avg1, avg2, 50, c, 'filled')
hold on
% set(gca,'yscale','log')
% set(gca,'xscale','log')
X = xlim;
Y = ylim;
plot([0 X(2)], [0 Y(2)], '--k')

% figure
% histogram(avg1(units), 12)
% hold on
% histogram(avg2(units), 12)

figure
plot([1 2], [avg1 avg2])
hold on
plot([1 2], [mean(avg1) mean(avg2)], 'k', 'LineWidth', 4)
set(gca,'yscale','log')
xlim([0.5 2.5])

units = int & su & grpidx;

figure
subplot(2, 1, 1)
stdshade(log10(fr.strd(units, :)), 0.3, 'k', fr.tstamps / 60 / 60, 3)
hold on
plot([datInfo.nsamps(3) datInfo.nsamps(3)] / 20000 / 60 / 60, ylim, '--k')
axis tight




nsamps = cumsum(datInfo.blockduration * fs);
nsamps = cumsum(datInfo.nsamps);

% strd firing rate for pyr and int
units = find(int & su & grpidx);
figure
data = fr.norm;
subplot(2, 1, 1)
plot(fr.tstamps / 60, (data(units, :))')
hold on
ydata = median((data(units, :)), 1);
plot(fr.tstamps / 60, ydata, 'k', 'LineWidth', 5)
% stdshade(data(units, :), 0.3, 'k', fr.tstamps / 60)
for i = 1 : length(nsamps) - 1
    plot([nsamps(i) nsamps(i)] / fs / 60, ylim, '--k')
end
axis tight
% ylim([0 5])
Y = ylim;
fill([states{2} fliplr(states{2})]' / 60, [Y(1) Y(1) Y(2) Y(2)],...
    'b', 'FaceAlpha', 0.15,  'EdgeAlpha', 0);
% title('PYR')
ylabel('norm MFR')

subplot(2, 1, 2)
data = fr.strd;
% units = find(su & grpidx);
plot(fr.tstamps / 60, (data(units, :))')
hold on
ydata = mean((data(units, :)), 1);
stdshade(data(units, :), 0.3, 'k', fr.tstamps / 60)
for i = 1 : length(nsamps) - 1
    plot([nsamps(i) nsamps(i)] / fs / 60, ylim, '--k')
end
axis tight
% set(gca, 'yscale', 'log')
Y = ylim;
fill([states{2} fliplr(states{2})]' / 60, [Y(1) Y(1) Y(2) Y(2)],...
    'b', 'FaceAlpha', 0.15,  'EdgeAlpha', 0);
% title('INT')
ylabel('MFR [Hz]')











subplot(2, 1, 2)
units = int & su & grpidx;
plot(fr.tstamps / 60 / 60, (fr.strd(units, :))')
hold on
ydata = mean((fr.strd(units, :)), 1);
plot(fr.tstamps / 60 / 60, ydata, 'k', 'LineWidth', 4)
set(gca, 'yscale', 'log')
plot([datInfo.nsamps(1) datInfo.nsamps(1)] / 20000 / 60 / 60, ylim, '--k')
axis tight


% mean + std of normalized fr
units = pyr & su & grpidx;
figure
stdshade(fr.strd(units, :), 0.3, 'k', fr.tstamps / 60, 3)
axis tight
hold on
nsamps = cumsum(datInfo.nsamps);
for i = 1 : length(nsamps) - 1
    plot([nsamps(i) nsamps(i)] / fs / 60, ylim, '--k')
end