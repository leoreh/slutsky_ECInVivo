function TL_catylzer_v2(cat , recDay)

%% FUNCTION: plots ripple analysis over recording time within user-inputted time range, using
% data concatenated across days

% INPUTS:
% [cat] : ripple data concatenated across days, outputted from TL_rippleWrapper_v2 -> TL_catMouse

% [recRange] : 2 element vector, each value indicating recording number for
% the start and stop of the analysis range (eg [1 2])
% [timeRange] : times for the onset and offsets of the analysis range, each
% pertaining to the recording time in days

% eg, (cat , [1 2] , [10.5 18]) will analysis data between 10:30AM on Day 1
% and 6PM on Day 2)

%% viewer and analyzer of data outputted from TL_rippleWrapper

% GUI that lets you scroll thru each ripple and possibly flag them

% Then option to open analysis pane
clearvars -global

%% Make state variables

% [stateList] : List of states that will refer to each corresponding value in states (eg.
% 1 = activeWake, 2 = quietWake...)
stateList = {'activeWake', 'quietWake' , 'lightSleep' , 'NREM' , 'REM' , 'N/REM' , 'bin' , 'unknown'};

% [statecols] : Color code for each state in stateList
stateCols = [0.8 0 0; 1 0.25 0.25; 0.25 0.25 0.8; 0 0 0.8; 0.25 0.8 0.25; 0.4 0.4 0.4; 0.7 0.7 0.7; 1 1 1];

%% Make Figure Panel and more variables

startDay = cat.start(recDay);
endDay = cat.stop(recDay);

xDistMin = 24*60*(endDay - startDay);

if xDistMin <= 720
    tik = 60;
end
if xDistMin > 720 & xDistMin < 1440
    tik = 120;
end
if xDistMin >= 1440
    tik = 240;
end
xaxTickz = [0 : tik : xDistMin - 5/60];
if xaxTickz(end) ~= xDistMin
    xaxTickz(end+1) = xDistMin;
end

% Make clever way to get ticklabels as time in day
% xlabelz = allHours

%% Extract ripples, mua, sua, states within user-inputted range

% Extract multi-unit activity
indx = cat.spt >= startDay & cat.spt < endDay;
currentdata.mu.spt = cat.spt(indx);
clear indx;

% Extract vals binned by second
indx = cat.dayT >= startDay & cat.dayT < endDay;
currentdata.ripBinnedHz = cat.ripHz(indx);
currentdata.muaBinnedHz = cat.MuaHz(indx);
currentdata.muaNonRipBinnedHz = cat.MuaNonRipHz(indx);
currentdata.sleep = cat.sleep(indx);
currentdata.dayT = cat.dayT(indx);

% % Extract single-unit activity
% % need to include nans in the state due to gaps between recordings...
% spt2 = cellfun(@(z) z(z>=startDay & z<endDay) - startDay , cat.su.spt , 'uniformoutput' , false);
% currentdata.su.spt = spt2(cell2mat(cellfun(@(z) ~isempty(z) , spt2 , 'uniformoutput' , false)));

% Extract ripplez
keep = cat.catRip.DayT(:,1) >= startDay & cat.catRip.DayT(:,2) < endDay;
fn = fieldnames(cat.catRip);
for f = 1 : length(fn)
currentdata.ripples.(fn{f})= cat.catRip.(fn{f})(keep,:);
end

%% Make Figure

plotWidth = 6 * xDistMin/(24*60);
F = figure('units' , 'inches' , 'position' , [0.5 , 0.5 , plotWidth + 2 ,  8] , ...
    'Name' , 'Control Panel' ,'Tag' , 'control_panel' , 'visible' , 'off');
set(F , 'renderer' , 'opengl');

%% Add ripple evoked spike rates over time

% Keep variables that are binned by xbin for possible use later on..

ax.pop.ripEvokedHz = axes('parent' , F , 'units' , 'inches' , 'position' , [1 6.5 plotWidth 1] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8);
hold on;
% set x axis limit
ax.pop.ripEvokedHz.XLim = [xaxTickz(1) , xaxTickz(end)];
ax.pop.ripEvokedHz.XTick = xaxTickz;
ax.pop.ripEvokedHzXLabel.String = 'Minutes';
ax.pop.ripEvokedHz.YLabel.String = 'Multi-Unit Firing Rate (Hz)';
ax.pop.ripEvokedHz.YAxis.Exponent = 0;

d = 24*60*60*1000*(currentdata.ripples.DayT(:,2) - currentdata.ripples.DayT(:,1));
currentdata.ripples.durationMs = d;


clear m p d;

% Calculate mean at each 5 minute interval during ripples and nonripple
% epochs
xbins = [0:60*5:ax.pop.ripEvokedHz.XLim(2)*60]; % 1 minute averages
secT = 24*60*60*(currentdata.dayT - startDay);

normRipSecT = 24*60*60*(currentdata.ripples.DayT - startDay);
% Calculate amount of time during each 5 minute bin that is during a ripple
carryOver = 0;
binned_ripTotSec = [];
for x = 1 : length(xbins) - 1
    totT = 0;
    i  = normRipSecT(:,1) >= xbins(x) & normRipSecT(:,2) < xbins(x+1);
    totT = totT + sum(normRipSecT(i,2) - normRipSecT(i,1));
    clear i;
    if x < length(xbins) - 1
        i = normRipSecT(:,1) >= xbins(x) & normRipSecT(:,1) < xbins(x+1) ...
            & normRipSecT(:,2) > xbins(x+1) & normRipSecT(:,2) < xbins(x+2);
        totT = totT + sum(xbins(x+1) - normRipSecT(i,1)); % should be 1 ripple max..
        carryOver = sum(normRipSecT(i,2) - xbins(x + 1));
        clear i;
    else
        i = normRipSecT(:,1) >= xbins(x) & normRipSecT(:,1) < xbins(x+1) ...
            & normRipSecT(:,2) > xbins(x+1) & normRipSecT(:,2) < xbins(x+1) + 300;
        totT = totT + sum(xbins(x+1) - normRipSecT(i,1)); % should be 1 ripple max..
        clear i;
    end
    totT = totT + carryOver;
    binned_ripTotSec(x) = totT;
end

% Calculate ripple and nonripple firing rates within each bin
for x = 1 : length(xbins) - 1
    
    i = secT >= xbins(x) & secT < xbins(x+1);
    binned_ripSpkCnt(x) = sum(currentdata.ripBinnedHz(i));
    
    % Need to calculate total time in this in in ripple to get correct Hz
    
    binned_ripHz(x) = binned_ripSpkCnt(x)/binned_ripTotSec(x);
    if binned_ripHz(x) > 0
    binned_ripHzSE(x) = std(currentdata.ripBinnedHz(i))/sqrt(sum(i));
    else
binned_ripHzSE(x) = 0;
    end
    binned_muSpkCnt(x) = sum(currentdata.muaBinnedHz(i));
    binned_nonRipHz(x) = sum(currentdata.muaNonRipBinnedHz(i))/(xbins(x+1)-xbins(x));

end

% if there are no ripple in an xbin, binned_ripHzSE will be NaN and this will
% disrupt the patch, so need to check if there are nans and then patch it
% in..patches (groups).
if sum(isnan(binned_ripHzSE) > 0)% take 5 min averages for plotting

    f = find(isnan(binned_ripHzSE));
    for ff = 1 : length(f)
        if f(ff) ~= 1
            if ff == 1 & f(1) ~= 1
                i = 1;
            end
            if ff == 1 & f(1) == 1
                i = 2;
            end
            if ff > 1
                i = f(ff-1)+1;
            end
            binI = i:f(ff)-1;
            patch(ax.pop.ripEvokedHz , 2.5 + [xbins(binI) , flip(xbins(binI))]/60 , [binned_ripHz(binI) - binned_ripHzSE(binI) , flip(binned_ripHz(binI) + binned_ripHzSE(binI))] , ...
                [0.4 0.4 0.4] , 'edgecolor' , 'none' , 'facealpha' , 0.5);
        end
    end
    if f(end) ~= length(binned_ripHzSE)
        binI = [f(end)+1:length(binned_ripHzSE)];
        patch(ax.pop.ripEvokedHz , 2.5 + [xbins(binI) , flip(xbins(binI))]/60 , [binned_ripHz(binI) - binned_ripHzSE(binI) , flip(binned_ripHz(binI) + binned_ripHzSE(binI))] , ...
            [0.4 0.4 0.4] , 'edgecolor' , 'none' , 'facealpha' , 0.5);
    end
else
    patch(ax.pop.ripEvokedHz , 2.5 + [xbins(1:end-1) , flip(xbins(1:end-1))]/60 , [binned_ripHz - binned_ripHzSE, flip(binned_ripHz + binned_ripHzSE)] , ...
        [0.4 0.4 0.4] , 'edgecolor' , 'none' , 'facealpha' , 0.5);
end
plot(ax.pop.ripEvokedHz , 2.5 + xbins(1:end-1)/60 , binned_ripHz , 'linestyle' , '-' , 'linewidth' , 0.5 , 'color' , 'c' , 'marker' , 'none');
patch(ax.pop.ripEvokedHz , [2.5 + xbins(1:end-1)/60 , 2.5 + xbins(end-1)/60 , 0] , [binned_nonRipHz , 0 , 0] , [0.4 0.4 0.4] , 'edgecolor' , 'none' , 'facealpha' , 0.5);
plot(ax.pop.ripEvokedHz , 2.5 + xbins(1:end-1)/60 , binned_nonRipHz , 'linestyle' , '-' , 'linewidth' , 0.5 , 'color' , 'k' , 'marker' , 'none');

val = find(cell2mat(cellfun(@(z) strcmp(z , 'bin') | strcmp(z , 'unknown') , stateList , 'uniformoutput' , false)));
exclMin = 24*60*60*(currentdata.dayT(ismember(currentdata.sleep , val)) - startDay);
exclbin = []
for x = 1 : length(xbins) - 1
    exclbin(x) = sum(exclMin >= xbins(x) & exclMin <= xbins(x+1)) > 0;
end
exclbin(70:110) = 1;
exclbin = logical(exclbin);

ymax = 1.25 * max([binned_ripHz(~exclbin) + binned_ripHzSE(~exclbin) , binned_nonRipHz(~exclbin)]);
ax.pop.ripEvokedHz.YLim = [0 ymax];

clear ymax riphzSE;

axbutt = ax.pop.ripEvokedHz.Position(2);

%% add patch of states directly beneath this plot

ax.pop.states = axes('parent' , F , 'units' , 'inches' , 'position' , [1 axbutt - 0.20 plotWidth 0.175] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'none' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8);
hold on;
ax.pop.states.XLim = [xaxTickz(1) , xaxTickz(end)];
ax.pop.states.YLim = [0 1];
tic
% Plot states per each second (True line width in inches is 'linewidth'/72
sval = currentdata.sleep(1:xbins(end));
ssec = 24*60*60*(currentdata.dayT(1:xbins(end)) - startDay);
for s = 1 : length(sval)
    if ~isnan(sval(s))
        t1 = ssec(s);
        if sval(s) <= 7
        patch(ax.pop.states , [t1-1,t1,t1,t1-1]/60 , [0 0 1 1] , stateCols(sval(s),:) , ...
            'edgecolor' , 'none' , 'facealpha' , 0.5);
        end
    end
end
toc

axbutt = ax.pop.states.Position(2);

%% Add patch for lights on/off periods
%kludge
startHr = 24*(startDay - floor(startDay));
endHr = 24*(endDay - floor(endDay));
endHr = endHr + 24*(floor(endDay) - floor(startDay))

if startHr < 8
    offstarts = [0 , (8-startHr):12:endHr-startHr];
end
if startHr >= 20
    offstarts = [0 , (24-startHr + 8):12:endHr-startHr];
end
if startHr >= 8 & startHr < 20
    offstarts = [20 - startHr : 12 : endHr - startHr + 20 - startHr];
end
 
ax.pop.lights = axes('parent' , F , 'units' , 'inches' , 'position' , [1 axbutt - 0.20 plotWidth 0.175] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'none' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8);
hold on;
ax.pop.lights.XLim = [xaxTickz(1) , xaxTickz(end)];
ax.pop.lights.YLim = [0 1];

cols = [0.3 0.3 0.3; 1 1 0.75];
for s = 1 : length(offstarts) - 1
    c = cols(~rem(s,2)+1,:);

        patch(ax.pop.lights , [offstarts(s) , offstarts(s+1) , offstarts(s+1) , offstarts(s)]*60 , [0 0 1 1] , c , ...
            'edgecolor' , 'none');
end

if offstarts(1) ~= 0
    patch(ax.pop.lights , [0 , offstarts(s) , offstarts(s) , 0]*60 , [0 0 1 1] , [1 1 0.75] , ...
            'edgecolor' , 'none');
end

axbutt = ax.pop.lights.Position(2);

%% Plot normalized firing from ripples and nonripple epochs
ax.pop.normFiring = axes('parent' , F , 'units' , 'inches' , 'position' , [1 axbutt - 0.75  plotWidth 0.5] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8 , 'XTickLabel' , []);
hold on;
% set x axis limit
ax.pop.normFiring.XLim = [xaxTickz(1) , xaxTickz(end)];
ax.pop.normFiring.XTick = xaxTickz;
% ax.pop.normFiring.XLabel.String = 'Minutes';
ax.pop.normFiring.YLabel.String = {'Normalized';'Firing'};
ax.pop.normFiring.YLim = [0 , 1];

% Plot ripple and nonripple firing rates scaled from 0 to 1 on the same
% axes
normRip = (binned_ripHz - min(binned_ripHz(~exclbin)))/(max(binned_ripHz(~exclbin)) - min(binned_ripHz(~exclbin)));
normnonRip = (binned_nonRipHz - min(binned_nonRipHz(~exclbin)))/(max(binned_nonRipHz(~exclbin)) - min(binned_nonRipHz(~exclbin)));
normRipCnt = (binned_ripSpkCnt - min(binned_ripSpkCnt(~exclbin)))/(max(binned_ripSpkCnt(~exclbin)) - min(binned_ripSpkCnt(~exclbin)));
plot(ax.pop.normFiring ,  2.5 + xbins(1:end-1)/60 , normnonRip , 'linestyle' , '-' , 'linewidth' , 0.5 , 'color' , 'k' , 'marker' , 'none');
plot(ax.pop.normFiring ,  2.5 + xbins(1:end-1)/60 , normRipCnt , 'linestyle' , '-' , 'linewidth' , 0.5 , 'color' , 'g' , 'marker' , 'none');
plot(ax.pop.normFiring ,  2.5 + xbins(1:end-1)/60 , normRip , 'linestyle' , '-' , 'linewidth' , 0.5 , 'color' , 'c' , 'marker' , 'none');


% Plot middle line..
plot(ax.pop.normFiring , ax.pop.normFiring.XLim , [0.5 0.5] , 'color' , [0.5 0.5 0.5] , 'linestyle' , ':' , 'linewidth' , 0.5);

axbutt = ax.pop.normFiring.Position(2);
%% Fraction of spikes during ripples

ax.pop.spikeFrac = axes('parent' , F , 'units' , 'inches' , 'position' , [1 axbutt - 0.75 plotWidth 0.5] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8 , 'XTickLabel' , []);
hold on;
ax.pop.spikeFrac.XLim = [xaxTickz(1) , xaxTickz(end)];
ax.pop.spikeFrac.XTick = xaxTickz;
% ax.pop.normFiring.XLabel.String = 'Minutes';
ax.pop.spikeFrac.YLabel.String = {'% MU Spikes';'in Ripple'};

ripSpkPerc = 100 * binned_ripSpkCnt ./ binned_muSpkCnt;

bar(ax.pop.spikeFrac , 2.5 + xbins(1:end-1)/60 , ripSpkPerc , 1  ,'facecolor' ,'k' , 'edgecolor' ,'none');

% add chance level (fraction of time in bin that is during  aripple)
tot = diff(xbins);
chance = 100 * binned_ripTotSec ./ tot;
plot(ax.pop.spikeFrac , 2.5 + xbins(1:end-1)/60 , chance , 'marker' , 'none' , 'linestyle' , '-' , 'color' , 'g');



axbutt = ax.pop.spikeFrac.Position(2);
%% Plot ripple rate over time

ax.pop.ripRate = axes('parent' , F , 'units' , 'inches' , 'position' , [1 axbutt - 0.75 plotWidth 0.5] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8);
hold on;

ax.pop.ripRate.XLim = [xaxTickz(1) , xaxTickz(end)];
ax.pop.ripRate.XTick = xaxTickz;
% ax.pop.ripRate.XLabel.String = 'Minutes';
ax.pop.ripRate.YLabel.String = {'Ripples';'/ Sec'};

p = 24*60*60*(currentdata.ripples.DayT(:,1) - startDay);
for x = 1 : length(xbins) - 1
    ripRate(x) = sum(p >= xbins(x) & p < xbins(x + 1))/300;
end
clear p;

bar(ax.pop.ripRate , 2.5 + xbins(1:end-1)/60 , ripRate , 1  ,'facecolor' ,'k' , 'edgecolor' ,'none');
% plot(ax.pop.normFiring ,  2.5 + xbins(1:end-1)/60 , normRip , 'linestyle' , '-' , 'linewidth' , 0.5 , 'color' , 'c' , 'marker' , 'none');

ripRate = ripRate;
clear ripRate;

axbutt = ax.pop.ripRate.Position(2);

%% Inter-ripple Interval

ax.pop.iri = axes('parent' , F , 'units' , 'inches' , 'position' , [1 axbutt - 0.75 plotWidth 0.5] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8);
hold on;
%

% With many ripples, probably best to do fraction of ripples occurring in a
% burst and then have another plot dots for each iri




ax.pop.iri.XLim = [xaxTickz(1) , xaxTickz(end)];
% ax.pop.iri.XTick = [0:60:ax.pop.ripEvokedHz.XLim(2)];
% ax.pop.iri.XLabel.String = 'Minutes';
ax.pop.iri.YLabel.String = 'IRI (ms)';
ax.pop.iri.YScale = 'linear';
ax.pop.iri.YDir = 'reverse';
ax.pop.iri.YLim = [0 750];
ax.pop.iri.YTick = [50,100,250 500];
ax.pop.iri.YTickLabel = {'25','100','250','500',[]};
%
p1 = 24*60*60*(currentdata.ripples.DayT(1:end-1,2)-startDay);
p2 = 24*60*60*(currentdata.ripples.DayT(2:end,1)-startDay);
iri = 1000 * (p2 - p1);
iriI = iri <= 500;
iriShow = iri(iriI);
x = p1(iriI);

% plot(ax.pop.iri , ax.pop.iri.XLim , repmat(cat.expInfo{1}.IRI,1,2) , 'linestyle' , ':' , 'linewidth' , 0.5 , 'color' , 'r' , 'marker' , 'none');
% % patch(ax.pop.iri , [0 ; p ; flip(p); 0] , [1000; iri ; flip(iri); 1000] , [0.75 0.75 0.75] , 'edgecolor' , 'none');
scatter(ax.pop.iri , x/60 , iriShow , 3 , 'k' , 'markerfacealpha' , 0.75 , 'markerfacecolor' , 'k' , 'markeredgecolor' , 'none' );% 'linestyle' , 'none' , 'marker' , 'o' , 'markersize' , 2 , 'markerfacecolor' , 'k'
% plot(ax.pop.iri , repmat(x',2,1)/60 , [1000*ones(1,length(x)); iriShow'] , 'linestyle' , '-' , 'color' , 'k' , 'marker' , 'none' , 'linewidth' , 0.15);

% plot yellow line at 250 ms just for reference for the 20percent line (so
% its easier to see the trend
plot(ax.pop.iri , ax.pop.iri.XLim , [250 250] , 'color' , 'y' , 'linestyle' , '-' , 'marker' , 'none' , 'linewidth' , 1);

p1v2 = p1(p1 < xbins(end));
[~ , ~ , iriPerBin] = histcounts(p1v2 , xbins);
u = [1:length(xbins)-1];
perc20 = [];
perc30 = [];
perc40 = [];
for uu = 1 : length(u)
    perc = [1:sum(iriPerBin == u(uu))]/sum(iriPerBin == u(uu));
    val = sort(iri(iriPerBin == u(uu)),'ascend');
    if length(val) > 20
    [~ , ind] = min(abs(perc - 0.20));
    perc20(uu) = val(ind);
    [~ , ind] = min(abs(perc - 0.30));
    perc30(uu) = val(ind);
    [~ , ind] = min(abs(perc - 0.40));
    perc40(uu) = val(ind);
    else
        perc20(uu) = nan;
        perc30(uu) = nan;
        perc40(uu) = nan;
    end
end

plot(ax.pop.iri , [xbins(1:end-1) - (xbins(2) - xbins(1))/2]/60 , perc20 , 'color' , 'c' , 'linestyle' , '-' , 'linewidth' , 1 , 'marker' , 'none');
% plot(ax.pop.iri , [xbins(1:end-1) - (xbins(2) - xbins(1))/2]/60 , perc30 , 'color' , 'g' , 'linestyle' , '-' , 'linewidth' , 1. , 'marker' , 'none');
% Overlay with median
axbutt = ax.pop.iri.Position(2);

%% Plot ripple duration over time
% 
% ax.pop.ripDur = axes('parent' , F , 'units' , 'inches' , 'position' , [1 axbutt - 0.75 plotWidth 0.5] , 'box' , 'off' , ...
%     'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8);
% hold on;
% 
% d = currentdata.ripples.durationMs;
% p = 24*60*60*(currentdata.ripples.DayT(:,1) - startDay);
% 
% for x = 1 : length(xbins) - 1
%     i = p >= xbins(x) & p < xbins(x+1);
%     dur(x) = mean(d(i));
%     durSE(x) = std(d(i))/sqrt(sum(i));
%     clear i;
% end
% clear d p;
% ax.pop.ripDur.XLim = [xaxTickz(1) , xaxTickz(end)];
% % ax.pop.ripDur.XTick = xaxTickz;
% ax.pop.ripDur.YLabel.String = {'Ripple';'Dur (ms)'};
% 
% bar(ax.pop.ripDur , 2.5 + xbins(1:end-1)/60 , dur , 1 , 'facecolor' , 'k' , 'edgecolor' , 'none' , 'linestyle' , 'none');
% 
% dur = dur;
% durSE = durSE;
% clear dur durSE;

F.Visible = 'on';

end
%% Fraction of single units participating in ripple
% 
% lI = currentdata.ripples.rangeSec;
% d = currentdata.ripples.durationMs;
% p = currentdata.ripples.peakSec;
% 
% for f = 1 : size(lI,1)
%     numU1(f) = sum(cell2mat(cellfun(@(z) sum(z >= lI(f,1) & z < lI(f,2) + d(f)/1000) >= 1 , su , 'uniformoutput' , false)));
%     numU2(f) = sum(cell2mat(cellfun(@(z) sum(z >= lI(f,1) & z < lI(f,2) + d(f)/1000) >= 2 , su , 'uniformoutput' , false)));
% end
% fracU1 = numU1/length(su);
% fracU2 = numU2/length(su);
% clear lI d spt numU1 numU2;
% 
% % Now average these two 5 minute bins, calculate averages and se
% for x = 1 : length(xbins) - 1
%     i = p >= xbins(x) & p < xbins(x+1);
%     m1(x) = mean(fracU1(i));
%     m2(x) = mean(fracU2(i));
%     se1(x) = std(fracU1(i))/sqrt(sum(i));
%     se2(x) = std(fracU2(i))/sqrt(sum(i));
%     clear i;
% end
% 
% ax.pop.fracU = axes('parent' , F , 'units' , 'inches' , 'position' , [1 ax.pop.ripEvokedHz.Position(2) - 4.5 4 0.5] , 'box' , 'off' , ...
%     'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'k' , 'fontname' , 'arial' , 'fontsize' , 8);
% hold on;
% 
% ax.pop.fracU.XLim = [xaxTickz(1) , xaxTickz(end)];
% ax.pop.fracU.XTick = xaxTickz
% ax.pop.fracU.YLabel.String = {'Fraction';'of Units'};
% 
% bar(ax.pop.fracU , 2.5 + xbins(1:end-1)/60 , m1 , 1 , 'facecolor' , 'k' , 'edgecolor' , 'none' , 'linestyle' , 'none');
% 
% bar(ax.pop.fracU , 2.5 + xbins(1:end-1)/60 , m2 , 1 , 'facecolor' , [0.75 0.75 0.75] , 'edgecolor' , 'none' , 'linestyle' , 'none');
% 
% su1spkFrac = m1;
% su1spkFracSE = se1;
% su2spkFrac = m2;
% su2spkFracSE = se2;
% 
% clear fracU1 fracU2 m1 m2 se1 se2 i
% %% Better than rastor plot may be for single units, fraction of spikes occurring during ripples with chance level marked as well
% 
% lI = currentdata.ripples.rangeSec;
% d = currentdata.ripples.durationMs;
% p = currentdata.ripples.peakSec;
% 
% ripSpkCnt = zeros(length(xbins) - 1 , length(su));
% totSpkCnt = zeros(size(ripSpkCnt));
% for x = 1 : length(xbins) - 1
%     
%     totSpkCnt(x,:) = cell2mat(cellfun(@(z) sum(z >= xbins(x) & z < xbins(x+1)) , su , 'uniformoutput' , false));
%     
%     i = p >= xbins(x) & p < xbins(x+1);
%     templI = lI(i,:);
%     for ii = 1 : size(templI , 1)
%         for s = 1 : length(su)
%             temp = sum(su{s} >= templI(ii,1) & su{s} < templI(ii,2));
%             ripSpkCnt(x,s) = ripSpkCnt(x,s) + temp;
%         end
%     end
% end
% 
% suRipFrac = ripSpkCnt ./ totSpkCnt;
% suRipFrac(isnan(suRipFrac)) = 0;
% 
% % Only plot cells that go above 0.05% at some point..
% temp = suRipFrac > 0.05;
% temp = sum(temp,1);
% i = temp == 0;
% suRipFrac(:,i) = nan;
% ax.pop.suRipFrac = axes('parent' , F , 'units' , 'inches' , 'position' , [1 ax.pop.ripEvokedHz.Position(2) - 5.25 4 0.5] , 'box' , 'off' , ...
%     'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8);
% hold on;
% ax.pop.suRipFrac.XLim = [xaxTickz(1) , xaxTickz(end)];
% % ax.pop.suRipFrac.XTick = [0:60:ax.pop.ripEvokedHz.XLim(2)];
% ax.pop.suRipFrac.YLim = [0 0.15];
% ax.pop.suRipFrac.YLabel.String = {'% Spikes';'in Ripple'};
% for s = 1 : size(suRipFrac , 2)
%     plot(ax.pop.suRipFrac , 2.5 + xbins(1:end-1)/60 , suRipFrac(: , s) , 'color' , rcols(s,:) , 'linestyle' , '-' , 'linewidth' , 0.5 , 'marker' , 'none');
% end
% 
% %% Add another patch of states below for clarity
% 
% ax.pop.states2 = axes('parent' , F , 'units' , 'inches' , 'position' , [1 ax.pop.ripEvokedHz.Position(2) - 5.45 4 0.175] , 'box' , 'off' , ...
%     'color' , 'none' , 'YColor' , 'none' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8);
% hold on;
% ax.pop.states2.XLim = [0 , ripples.expInfo.stopHr*60];
% % ax.pop.states2.XTick = [0 : 60 : ripples.expInfo.stopHr*60];
% ax.pop.states2.YLim = [0 1];
% tic
% % patch states for each second...
% for s = 1 : length(states)
%     val = states(s);
%     if ~isnan(val)
%         patch(ax.pop.states2 , [s-1,s,s,s-1]/60 , [0 0 1 1] , stateCols(states(s),:) , ...
%             'edgecolor' , 'none' , 'facealpha' , 0.5);
%     end
% end
% toc
% 
% %% Add x axis at the bottom
% 
% ax.pop.xax = axes('parent' , F , 'units' , 'inches' , 'position' , [1 ax.pop.ripEvokedHz.Position(2) - 5.5 4 0.175] , 'box' , 'off' , ...
%     'color' , 'none' , 'YColor' , 'none' , 'XColor' , 'k' , 'fontname' , 'arial' , 'fontsize' , 8);
% hold on;
% ax.pop.xax.XLim = [xaxTickz(1), xaxTickz(end)];
% ax.pop.xax.XTick = xaxTickz;
% ax.pop.xax.XLabel.String = 'Minutes';
% 
% %% Analyses Plots
% 
% % "Ripple-Triggered Averages" with different event types
% 
% % (1) Ripple bursts, where the cell fires at least 2 spikes in at least 1
% % of the ripples
% 
% % Find Ripple Bursts and the start and end of them
% p = currentdata.ripples.peakSec(1:end-1);
% iri = 1000 * (currentdata.ripples.peakSec(2:end) - p);
% 
% a = find(iri<500);
% lI = currentdata.ripples.rangeSec;
% % Find start and end times of each ripple bursts
% for aa = 1 : length(a)
%     ripBurstLimSec(aa,:) = [lI(a(aa),1) , lI(a(aa)+1,2)];
% end
% 
% % Only keep the bursts that are totally in NREM
% 
% %% Make Figure visible
% set(F,'visible' , 'on');
% %% single unit rastor plot (v slow maybe just do ripples)
% % TL determined this is not such a useful plot..
% 
% % ax.suRipRastor = axes('parent' , control_panel , 'units' , 'inches' , 'position' , [1 ax.ripEvokedHz.Position(2) - 5.0 4 1] , 'box' , 'off' , ...
% %     'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'k' , 'fontname' , 'arial' , 'fontsize' , 8 , 'XLim' , [0 6*60] , 'XTick' , [0:60:6*60] , 'TickLength' , [0 0]);
% % hold on;
% % ax.suRipRastor.YLabel.String = 'Unit';
% % ax.suRipRastor.YLim = [0 length(spikes.times) - 1];
% % ax.suRipRastor.XLim = [0 , ripples.expInfo.stopHr*60];
% % ax.suRipRastor.XTick = [0:60:ax.ripEvokedHz.XLim(2)];
% %
% % lI = currentData.rangeSec;
% % for s = 1 : length(spikes.times)
% %     spt = spikes.times{s};
% %     ripEvoked = [];
% %     for i = 1 : size(lI,1)
% %         temp = spt(spt >= lI(i,1) & spt < lI(i,2));
% %         ripEvoked = [ripEvoked; temp];
% %         clear temp;
% %     end
% %     plot(ax.suRipRastor , repmat(ripEvoked'/60,2,1) , repmat([s-1,s]',1,length(ripEvoked)) , 'linestyle' , '-' , 'color' , rcols(s,:) , 'marker' , 'none' , 'linewidth' ,0.5);
% % end
% 
% %% Callback functionz and other nested functions
% 
% % function cb_states(~ , evt)
% % global stateList ripInclude data;
% % type = evt.Source.Tag;
% % val = evt.Source.Value;
% %
% % state = find(strcmp(stateList , type));
% % ripInclude(data.ripList.states == state) = val;
% % clear type val state;
% 
% % function cb_rescale(~ , ~)
% %
% % function cb_sort(~ , evt)
% % global h data ripInclude;
% % temp = fieldnames(h.sort);
% % % set all colors to grey and values to zero
% % for t = 1 : length(temp)
% %     h.sort.(temp{t}).BackgroundColor = [0.94 0.94 0.94];
% %     h.sort.(temp{t}).UserData = 0;
% % end
% % clear temp t;
% % h.sort.(evt.Source.Tag).BackgroundColor = [0.94 0.94 0.5];
% % h.sort.(evt.Source.Tag).UserData = 1;
% %
% % function cb_ripplyzeIt(~ , evt)
% % global h data ripInclude stateList control_panel;
% %
% % for n = 1 : 3
% %     h.(evt.Source.Tag).BackgroundColor = [0.5 0.94 0.94];
% %     pause(0.1);
% %     h.(evt.Source.Tag).BackgroundColor = [0.94 0.94 0.94];
% %     pause(0.1);
% % end
% % % restore logical indices to include all the data
% % ripInclude = logical(ones(length(data.ripList.peakPosition) , 1));
% %
% % % Only include the selected states
% % hfn = fieldnames(h.state);
% % for hh = 1 : length(hfn)
% %     type = h.state.(hfn{hh}).Tag;
% %     val = h.state.(hfn{hh}).Value;
% %     state = find(strcmp(stateList , type));
% %     ripInclude(data.ripList.states == state) = val;
% % end
% %
% % % Now reorder ripInclude and the actual data sets
% % % based on value of h.descend and user-selected category
% % if h.descend.Value == 1
% %     direction = 'descend';
% % else
% %     direction = 'ascend';
% % end
% % % Find which sorting variable has been chosen
% % % h.sort
% % temp = structfun(@(z) z.UserData == 1 , h.sort , 'uniformoutput' , false);
% % fn = fieldnames(temp);
% % categoryIndx = cell2mat(struct2cell(temp));
% % [~ , order] = sort(data.ripList.(fn{categoryIndx})(ripInclude) , direction);
% %
% % fn = fieldnames(data.ripList);
% % for f = 1 : length(data.fn)
% %     temp = data.ripList.(fn{f})(ripInclude);
% %     currentData.(fn{f}) = temp(order);
% % end
% % setappdata(control_panel , 'currentData' , currentData);
% %
% % h.index.String = '1';
% % % Clear all axes
% % makeRipplePlots(currentData , str2num(h.index.String));
% %
% % function cb_plus(~ , ~)
% % global h control_panel
% %
% % currentData = getappdata(control_panel , 'currentData');
% %
% % % Change h.index.Value (or string)
% % val = str2num(h.index.String);
% % fn = fieldnames(currentData);
% % if val ~= size(currentData.(fn{1}),1)
% %     h.index.String = num2str(str2num(h.index.String) + 1);
% % end
% %
% % % Run Plotting function
% % makeRipplePlots(currentData , str2num(h.index.String));
% %
% % function cb_minus(~ , ~)
% % global h control_panel
% %
% % currentData = getappdata(control_panel , 'currentData');
% %
% % % Change h.index.Value (or string)
% % val = str2num(h.index.String);
% % if val ~= 1
% %     h.index.String = num2str(val - 1);
% % end
% %
% % % Run Plotting function
% % makeRipplePlots(currentData , str2num(h.index.String));
% %
% % % function makeRipplePlots(d , indx)
% % % % clear all axes..
% % %
% % % %
% % % %    c = get(ax , 'children');
% % % %             for cc = 1 : length(c)
% % % %                 delete(c(cc));
% % % %             end
% %
% %
% 
% 
% function [limz] = findLims(cellArray)
% limz = [min(cell2mat(cellfun(@(z) min(z) , cellArray , 'uniformoutput' , false))) , max(cell2mat(cellfun(@(z) max(z) , cellArray , 'uniformoutput' , false)))];
% 
% % function [butt] = scaleAx(~ , evt)
% % global ax butt
% %
% % % left click
% % if evt.Button == 1
% %     butt = evt.IntersectionPoint(1);
% % end
% %
% % % right click
% % if evt.Button == 3
% %     if exist('butt') == 1
% %         if ~isempty(butt)
% %             if evt.IntersectionPoint(2) > butt(1)
% %                 butt(2) = evt.IntersectionPoint(2);
% %                 fn = fieldnames(ax.pop)
% %                 for f = 1 : length(fn)
% %                     d = round(diff(butt)/2 + butt(1),0);
% %                     ax.pop.(fn{f}).XLim = butt;
% %                 end
% %                 ax.pop.xax.XTick = [butt(1) , d , butt(2)];
% %             end
% %         end
% %     end
% %     butt = [];
% % end
% %
% function  cb_xaxTickz(~,~)
% global ax xaxTickz
% fn = fieldnames(ax.pop)
% for f = 1 : length(fn)
%     ax.pop.(fn{f}).XLim = [xaxTickz(1) , xaxTickz(end)];
% end
% ax.pop.xax.XTick = xaxTickz;


% 
% 
% 
