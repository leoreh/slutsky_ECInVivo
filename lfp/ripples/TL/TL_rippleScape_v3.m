function TL_rippleScape_v3(cat , recRange , timeRange)

%% FUNCTION: plots ripple analysis over recording time within user-inputted time range, using
% data concatenated across days

% INPUTS:
% [cat] : ripple data concatenated across days, outputted from TL_rippleWrapper_v2 -> TL_catMouse
% [recRange] : 2 element vector, each value indicating recording number for
% the start and stop of the analysis range (eg [1 2])
% [timeRange] : times for the onset and offsets of the analysis range, each
% pertaining to the time in recording indices recRange(1) and recRange(2).
% decimal points refer to fraction of the hour
gca
% eg, (cat , [1 2] , [10.5 18]) will analysis data between 10:30AM on Day 1
% and 6PM on Day 2)

%% viewer and analyzer of data outputted from TL_rippleWrapper

% GUI that lets you scroll thru each ripple and possibly flag them

% Then option to open analysis pane
clearvars -global

%% Make state variables

% [stateList] : List of states that will refer to each corresponding value in states (eg.
% 1 = activeWake, 2 = quietWake...)
stateList = {'activeWake', 'quietWake' , 'lightSleep' , 'NREM' , 'REM' , 'bin' , 'unknown'};

% [statecols] : Color code for each state in stateList
stateCols = [0.8 0 0; 1 0.25 0.25; 0.25 0.25 0.8; 0 0 0.8; 0.25 0.8 0.25; 0.4 0.4 0.4; 0.7 0.7 0.7];

%% Make Figure Panel and more variables


% % [rcols] : Color code for each single unit
j = repmat(jet,100,1);
rcols = j([1:length(cat.su.spt)],:);

% % [lfpSampRate] : samples per second in lfp recordings
% lfpSampRate = ripples.expInfo.lfpSampRate;

% % [bufferMs] : time in ms saved as a buffer before and after each ripple in
% % the ripples structure
% bufferMs = ripples.expInfo.bufferMs;

% [xDist] : total length in minutes of the analysis window
% [xaxTickz] : x axis tix and limits used for plotting. Data will be binned
% by 5 minutes, so subtract 5 min from the end
xDist = 60 * ((24 - timeRange(1)) + (diff(recRange)-1)*24 + timeRange(2));
if xDist <= 720
    tik = 60;
end
if xDist > 720 & xDist < 1440
    tik = 120;
end
if xDist >= 1440
    tik = 240;
end
xaxTickz = [0 : tik : xDist - 5/60];
if xaxTickz(end) ~= xDist
    xaxTickz(end+1) = xDist;
end

% Make clever way to get ticklabels as time in day
% xlabelz = allHours

%% Extract ripples, mua, sua, states within user-inputted range

startSec  = 60 * 60 * (24 - cat.expInfo{1}.RecStartTime + (recRange(1) - 1 - 1)*24 + timeRange(1));
stopSec = startSec + xDist*60;

% Extract multi-unit activity
indx = cat.mu.spt >= startSec & cat.mu.spt < stopSec;
currentdata.mu.spt = cat.mu.spt(indx) - startSec;
currentdata.mu.recNum = cat.mu.recNum(indx);
clear indx;

% Extract single-unit activity
% need to include nans in the state due to gaps between recordings...
spt2 = cellfun(@(z) z(z>=startSec & z<stopSec) - startSec , cat.su.spt , 'uniformoutput' , false);
indx = cell2mat(cellfun(@(z) ~isempty(z) , spt2 , 'uniformoutput' , false));
currentdata.su.spt = spt2(indx);
currentdata.su.recNum = cat.su.recNum(indx);
clear spt2 indx;

% Extract ripplez
keep = cat.ripples.rangeSec(:,1) >= startSec & cat.ripples.rangeSec(:,2) < stopSec;
fn = fieldnames(cat.ripples);
for f = 1 : length(fn)
currentdata.ripples.(fn{f})= cat.ripples.(fn{f})(keep,:);
end
currentdata.ripples.peakSec = currentdata.ripples.peakSec - startSec;
currentdata.ripples.rangeSec = currentdata.ripples.rangeSec - startSec;

temp = cat.gap - startSec;
keep = temp(:,1) > 0 & temp(:,1) < xDist*60 | temp(:,2) > 0 & temp(:,2) < xDist*60;
currentdata.gaps = temp(keep,:);

% Extract statez
indx = floor(cat.state.sec) >= startSec & cat.state.sec < stopSec;
currentdata.state.val = cat.state.val(indx);
currentdata.state.sec = cat.state.sec(indx) - startSec;
currentdata.state.val(currentdata.state.val == 0) = nan;

% indx of time bins coinciding with gaps;
xbins = [0:60*5:xaxTickz(end)*60]; % 1 minute averages
excl = [];
for g = 1 : size(currentdata.gaps , 1)
    for x = 1 : length(xbins) - 1
        excl(g,x) = currentdata.gaps(g,1) >= xbins(x) & currentdata.gaps(g,1) < xbins(x+1) ...
            | currentdata.gaps(g,2) >= xbins(x) & currentdata.gaps(g,2) < xbins(x+1) ...
            | currentdata.gaps(g,1) < xbins(x) & currentdata.gaps(g,2) > xbins(x+1);
    end
end
excl = sum(excl,1)>0;
%+/- 1 bin
excl(find(diff(excl) == 1) - 1) = 1;
excl(find(diff(excl) == -1) + 1) = 1;
%% Make Figure

plotWidth = 6 * xDist/(24*60);
F = figure('units' , 'inches' , 'position' , [0.5 , 0.5 , plotWidth + 2 ,  8] , 'PaperPosition' , [0, 0 , plotWidth + 2 ,  8] , ...
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

% Calculate multiunit firing rate during each ripple
m = currentdata.ripples.muSpkTimesMs;
d = currentdata.ripples.durationMs;
for r = 1 : size(currentdata.ripples.rangeSec , 1)
    hz = sum(m{r} >= 0 & m{r} < d(r)) ...
        / (d(r)/1000);
    currentdata.ripples.ripHz(r,1) = hz;
    clear hz;
end
clear m d p;

% Calculate mean at each 5 minute interval during ripples and nonripple
% bouts
xbins = [0:60*5:ax.pop.ripEvokedHz.XLim(2)*60]; % 1 minute averages
ripHz = currentdata.ripples.ripHz;
p = currentdata.ripples.peakSec;
d = currentdata.ripples.durationMs;
m = currentdata.mu.spt;
for x = 1 : length(xbins) - 1
    if x == length(xbins) - 1
    end
    i = p >= xbins(x) & ...
        p < xbins(x+1);
    ripHzMean(x) = mean(ripHz(i));
    ripHzSE(x) = std(ripHz(i))/sqrt(sum(i));
    sumI(x) = sum(i);
    % find nonripple evoked firing rate
    nonRipHz(x) = (sum(m >= xbins(x) & m < xbins(x+1)) ...
        -  (sum(ripHz(i)) * sum(d(i)/1000))) ...
        / 300;
end
clear p d;

% take 5 min averages for plotting
% if there are no ripple in an xbin, ripHzSE will be NaN and this will
% disrupt the patch, so need to check if there are nans and then patch it
% in..patches (groups).
if sum(isnan(ripHzSE) > 0)
    f = find(isnan(ripHzSE));
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
            patch(ax.pop.ripEvokedHz , 2.5 + [xbins(binI) , flip(xbins(binI))]/60 , [ripHzMean(binI) - ripHzSE(binI) , flip(ripHzMean(binI) + ripHzSE(binI))] , ...
                [0.4 0.4 0.4] , 'edgecolor' , 'none' , 'facealpha' , 0.5);
        end
    end
    if f(end) ~= length(ripHzSE)
        binI = [f(end)+1:length(ripHzSE)];
        patch(ax.pop.ripEvokedHz , 2.5 + [xbins(binI) , flip(xbins(binI))]/60 , [ripHzMean(binI) - ripHzSE(binI) , flip(ripHzMean(binI) + ripHzSE(binI))] , ...
            [0.4 0.4 0.4] , 'edgecolor' , 'none' , 'facealpha' , 0.5);
    end
else
    patch(ax.pop.ripEvokedHz , 2.5 + [xbins(1:end-1) , flip(xbins(1:end-1))]/60 , [ripHzMean - ripHzSE, flip(ripHzMean + ripHzSE)] , ...
        [0.4 0.4 0.4] , 'edgecolor' , 'none' , 'facealpha' , 0.5);
end
plot(ax.pop.ripEvokedHz , 2.5 + xbins(1:end-1)/60 , ripHzMean , 'linestyle' , '-' , 'linewidth' , 0.5 , 'color' , 'c' , 'marker' , 'none');
patch(ax.pop.ripEvokedHz , [2.5 + xbins(1:end-1)/60 , 2.5 + xbins(end-1)/60 , 0] , [nonRipHz , 0 , 0] , [0.4 0.4 0.4] , 'edgecolor' , 'none' , 'facealpha' , 0.5);
plot(ax.pop.ripEvokedHz , 2.5 + xbins(1:end-1)/60 , nonRipHz , 'linestyle' , '-' , 'linewidth' , 0.5 , 'color' , 'k' , 'marker' , 'none');

binD.ripHzMean = ripHzMean;
binD.ripHzSE = ripHzSE;
binD.nonRipHz = nonRipHz;
clear ripHzMean riphzSE nonRipHz;

axbutt = ax.pop.ripEvokedHz.Position(2);

%% add patch of states directly beneath this plot

ax.pop.states = axes('parent' , F , 'units' , 'inches' , 'position' , [1 axbutt - 0.20 plotWidth 0.175] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'none' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8);
hold on;
ax.pop.states.XLim = [xaxTickz(1) , xaxTickz(end)];
ax.pop.states.YLim = [0 1];
tic
% Plot states per each second (True line width in inches is 'linewidth'/72
sval = currentdata.state.val;
ssec = currentdata.state.sec;
for s = 1 : length(sval)
    if ~isnan(sval(s))
        t1 = ssec(s);
        patch(ax.pop.states , [t1-1,t1,t1,t1-1]/60 , [0 0 1 1] , stateCols(sval(s),:) , ...
            'edgecolor' , 'none' , 'facealpha' , 0.5);
    end
end
toc

axbutt = ax.pop.states.Position(2);

%% Add patch for lights on/off periods
%kludge
lightsOn = 7;
% lightsOn = cat.expInfo{1}.LightsTime;
switchz = 60*60*(lightsOn - timeRange(1):12:(stopSec-startSec)/60/60+12);
if switchz(end) < xaxTickz(end)
    switchz(end+1) = xaxTickz(end);
end
switchz = switchz(switchz>0);
if switchz(1) ~= 0
    switchz = [0 , switchz];
end

ax.pop.lights = axes('parent' , F , 'units' , 'inches' , 'position' , [1 axbutt - 0.20 plotWidth 0.175] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'none' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8);
hold on;
ax.pop.lights.XLim = [xaxTickz(1) , xaxTickz(end)];
ax.pop.lights.YLim = [0 1];

% Check if period begins during light or dark and specify colors accordingly
if cat.expInfo{1}.RecStartTime >= lightsOn & cat.expInfo{1}.RecStartTime < lightsOn + 12
cols = [1 1 0.75; 0.1 0.1 0.1];
else
    cols = [0.1 0.1 0.1; 1 1 0.75];
end

for s = 1 : length(switchz) - 1
    c = cols(~rem(s,2)+1,:);
    t1 = switchz(s);
    t2 = switchz(s+1);
        patch(ax.pop.lights , [t1,t2,t2,t1]/60 , [0 0 1 1] , c , ...
            'edgecolor' , 'none' , 'facealpha' , 0.5);
end

axbutt = ax.pop.lights.Position(2);
% 
%% Plot normalized firing from ripples and nonripple bouts
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
normRip = (binD.ripHzMean - min(binD.ripHzMean(~excl)))/(max(binD.ripHzMean(~excl)) - min(binD.ripHzMean(~excl)));
normnonRip = (binD.nonRipHz - min(binD.nonRipHz(~excl)))/(max(binD.nonRipHz(~excl)) - min(binD.nonRipHz(~excl)));
plot(ax.pop.normFiring ,  2.5 + xbins(1:end-1)/60 , normRip , 'linestyle' , '-' , 'linewidth' , 0.5 , 'color' , 'c' , 'marker' , 'none');
plot(ax.pop.normFiring ,  2.5 + xbins(1:end-1)/60 , normnonRip , 'linestyle' , '-' , 'linewidth' , 0.5 , 'color' , 'k' , 'marker' , 'none');

% Plot middle line..
p = plot(ax.pop.normFiring , ax.pop.normFiring.XLim , [0.5 0.5] , 'color' , [0.5 0.5 0.5] , 'linestyle' , ':' , 'linewidth' , 0.5);
uistack(p,'bottom');
p = plot(ax.pop.normFiring , ax.pop.normFiring.XLim , [0 0] , 'color' , [0.5 0.5 0.5] , 'linestyle' , ':' , 'linewidth' , 0.5);
uistack(p,'bottom');
p = plot(ax.pop.normFiring , ax.pop.normFiring.XLim , [1 1] , 'color' , [0.5 0.5 0.5] , 'linestyle' , ':' , 'linewidth' , 0.5);
uistack(p,'bottom');

axbutt = ax.pop.normFiring.Position(2);
%% Fraction of spikes during ripples

ax.pop.spikeFrac = axes('parent' , F , 'units' , 'inches' , 'position' , [1 axbutt - 0.75 plotWidth 0.5] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8 , 'XTickLabel' , []);
hold on;
ax.pop.spikeFrac.XLim = [xaxTickz(1) , xaxTickz(end)];
ax.pop.spikeFrac.XTick = xaxTickz;
% ax.pop.normFiring.XLabel.String = 'Minutes';
ax.pop.spikeFrac.YLabel.String = {'% MU Spikes';'in Ripple'};

lI = currentdata.ripples.rangeSec;
p = currentdata.ripples.peakSec;

ripSpkCnt = zeros(length(xbins) - 1 , 1);
totSpkCnt = zeros(size(ripSpkCnt));
carryover = 0;
for x = 1 : length(xbins) - 1
    
    totSpkCnt(x,1) = sum(m >= xbins(x) & m < xbins(x+1));
    ripSpkCnt(x,1) = carryover;
    carryover = 0;
    
    i = p >= xbins(x) & p < xbins(x+1);
    templI = lI(i,:);
    
    for ii = 1 : size(templI , 1)
        if ii == size(templI,1) && templI(2) > xbins(x+1)
            temp = sum(m >= templI(ii,1) & m < xbins(x+1));
            carryover = sum(m >= xbins(x+1) & m < templI(ii,2));
        else
            temp = sum(m >= templI(ii,1) & m < templI(ii,2));
        end
        ripSpkCnt(x,1) = ripSpkCnt(x,1) + temp;
    end
end

ripSpkPerc = 100 * ripSpkCnt ./ totSpkCnt;
bar(ax.pop.spikeFrac , 2.5 + xbins(1:end-1)/60 , ripSpkPerc , 1  ,'facecolor' ,'k' , 'edgecolor' ,'none');

binD.ripSpkCnt = ripSpkCnt;
binD.totSpkCnt = totSpkCnt;
binD.ripSpkPerc = ripSpkPerc;
clear ripSpkCnt totspkCnt ripSpkPerc;

axbutt = ax.pop.spikeFrac.Position(2);
%% Plot ripple rate over time

ax.pop.ripRate = axes('parent' , F , 'units' , 'inches' , 'position' , [1 axbutt - 0.75 plotWidth 0.5] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8);
hold on;

ax.pop.ripRate.XLim = [xaxTickz(1) , xaxTickz(end)];
ax.pop.ripRate.YLabel.String = {'Ripples';'/ Sec'};

p = currentdata.ripples.peakSec;
for x = 1 : length(xbins) - 1
    ripRate(x) = sum(p >= xbins(x) & p < xbins(x + 1))/300;
end
clear p;

bar(ax.pop.ripRate , 2.5 + xbins(1:end-1)/60 , ripRate , 1  ,'facecolor' ,'k' , 'edgecolor' ,'none');
% plot(ax.pop.normFiring ,  2.5 + xbins(1:end-1)/60 , normRip , 'linestyle' , '-' , 'linewidth' , 0.5 , 'color' , 'c' , 'marker' , 'none');

binD.ripRate = ripRate;
clear ripRate;

axbutt = ax.pop.ripRate.Position(2);

%% Plot Power

ax.pop.ripPower = axes('parent' , F , 'units' , 'inches' , 'position' , [1 axbutt - 0.75 plotWidth 0.5] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8);
hold on;
%
ax.pop.ripPower.XLim = [xaxTickz(1) , xaxTickz(end)];
ax.pop.ripPower.YLabel.String = 'Power (a.u.)';

p = currentdata.ripples.peakSec;
pnp = currentdata.ripples.peakNormalizedPower;
for x = 1 : length(xbins) - 1
    i = p >= xbins(x) & p < xbins(x + 1); 
    ripPower(x) = mean(pnp(i));
    clear i;
end
clear p pnp;

bar(ax.pop.ripPower , 2.5 + xbins(1:end-1)/60 , ripPower , 1  ,'facecolor' ,'k' , 'edgecolor' ,'none');

axbutt = ax.pop.ripPower.Position(2);
%% Inter-ripple Interval

ax.pop.iri = axes('parent' , F , 'units' , 'inches' , 'position' , [1 axbutt - 0.75 plotWidth 0.5] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8);
hold on;
%
ax.pop.iri.XLim = [xaxTickz(1) , xaxTickz(end)];
% ax.pop.iri.XTick = [0:60:ax.pop.ripEvokedHz.XLim(2)];
% ax.pop.iri.XLabel.String = 'Minutes';
ax.pop.iri.YLabel.String = 'IRI (ms)';
ax.pop.iri.YScale = 'linear';
ax.pop.iri.YDir = 'reverse';
ax.pop.iri.YLim = [0 750];
ax.pop.iri.YTick = [cat.expInfo{1}.IRI,250 500 750];
ax.pop.iri.YTickLabel = {num2str(cat.expInfo{1}.IRI),'250','500',[]};
%
p = currentdata.ripples.peakSec(1:end-1);
iri = 1000 * (currentdata.ripples.peakSec(2:end) - p);


plot(ax.pop.iri , ax.pop.iri.XLim , repmat(cat.expInfo{1}.IRI,1,2) , 'linestyle' , ':' , 'linewidth' , 0.5 , 'color' , 'r' , 'marker' , 'none');
% patch(ax.pop.iri , [0 ; p ; flip(p); 0] , [1000; iri ; flip(iri); 1000] , [0.75 0.75 0.75] , 'edgecolor' , 'none');
plot(ax.pop.iri , repmat(p'/60,2,1) , [1000*ones(1,length(p)); iri'] , 'linestyle' , '-' , 'color' , 'k' , 'marker' , 'none' , 'linewidth' , 0.15);

axbutt = ax.pop.iri.Position(2);

%% Plot ripple duration over time
% 
% ax.pop.ripDur = axes('parent' , F , 'units' , 'inches' , 'position' , [1 axbutt - 3.75 plotWidth 0.5] , 'box' , 'off' , ...
%     'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8);
% hold on;
% 
% d = currentdata.ripples.durationMs;
% p = currentdata.ripples.peakSec;
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
% binD.dur = dur;
% binD.durSE = durSE;
% clear dur durSE;

%% Fraction of single units participating in ripple

% Only include units recorded on the same day as each ripple
rS = currentdata.ripples.rangeSec;
d = currentdata.ripples.durationMs;
p = currentdata.ripples.peakSec;
rr = currentdata.ripples.recNum;
su = currentdata.su.spt;
sr = currentdata.su.recNum;

ur = unique(rr);
fracU1 = [];
fracU2 = [];
for u = 1 : length(ur)
    L = find(rr == ur(u));
    for el = 1 : length(L)
    fracU1(end+1) = sum(cell2mat(cellfun(@(z) sum(z >= rS(L(el),1) & z < rS(L(el),2)) >= 1 , [su(sr == ur(u))] , 'uniformoutput' , false)))/sum(sr == ur(u));
    fracU2(end+1) = sum(cell2mat(cellfun(@(z) sum(z >= rS(L(el),1) & z < rS(L(el),2)) >= 2 , [su(sr == ur(u))] , 'uniformoutput' , false)))/sum(sr == ur(u));
    end
end
% for f = 1 : size(lI,1)
%     numU1(f) = sum(cell2mat(cellfun(@(z) sum(z >= lI(f,1) & z < lI(f,2)) >= 1 , su , 'uniformoutput' , false)));
%     numU2(f) = sum(cell2mat(cellfun(@(z) sum(z >= lI(f,1) & z < lI(f,2)) >= 2 , su , 'uniformoutput' , false)));
% end
% end
% fracU1 = numU1/length(su);
% fracU2 = numU2/length(su);
% clear lI d spt numU1 numU2;

% Now average these two 5 minute bins, calculate averages and se
for x = 1 : length(xbins) - 1
    i = p >= xbins(x) & p < xbins(x+1);
    m1(x) = mean(fracU1(i));
    m2(x) = mean(fracU2(i));
    se1(x) = std(fracU1(i))/sqrt(sum(i));
    se2(x) = std(fracU2(i))/sqrt(sum(i));
    clear i;
end

ax.pop.fracU = axes('parent' , F , 'units' , 'inches' , 'position' , [1 axbutt - 0.75 plotWidth 0.5] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8);
hold on;

ax.pop.fracU.XLim = [xaxTickz(1) , xaxTickz(end)];
ax.pop.fracU.XTick = xaxTickz
ax.pop.fracU.YLabel.String = {'Fraction';'of Units'};

bar(ax.pop.fracU , 2.5 + xbins(1:end-1)/60 , m1 , 1 , 'facecolor' , 'k' , 'edgecolor' , 'none' , 'linestyle' , 'none');

bar(ax.pop.fracU , 2.5 + xbins(1:end-1)/60 , m2 , 1 , 'facecolor' , [0.75 0.75 0.75] , 'edgecolor' , 'none' , 'linestyle' , 'none');

binD.su1spkFrac = m1;
binD.su1spkFracSE = se1;
binD.su2spkFrac = m2;
binD.su2spkFracSE = se2;

clear fracU1 fracU2 m1 m2 se1 se2 i

axbutt = ax.pop.fracU.Position(2);
%% Better than rastor plot may be for single units, fraction of spikes occurring during ripples with chance level marked as well

rS = currentdata.ripples.rangeSec;
d = currentdata.ripples.durationMs;
p = currentdata.ripples.peakSec;
su = currentdata.su.spt;

ripSpkCnt = zeros(length(xbins) - 1 , length(su));
totSpkCnt = zeros(size(ripSpkCnt));
for x = 1 : length(xbins) - 1
    
    totSpkCnt(x,:) = cell2mat(cellfun(@(z) sum(z >= xbins(x) & z < xbins(x+1)) , su , 'uniformoutput' , false));
    
    i = p >= xbins(x) & p < xbins(x+1);
    temprS = rS(i,:);
    for ii = 1 : size(temprS , 1)
        for s = 1 : length(su)
            temp = sum(su{s} >= temprS(ii,1) & su{s} < temprS(ii,2));
            ripSpkCnt(x,s) = ripSpkCnt(x,s) + temp;
        end
    end
end

suRipPerc = 100 * ripSpkCnt ./ totSpkCnt;
suRipPerc(isnan(suRipPerc)) = 0;

% Only plot cells that go above 0.05% at some point..
temp = suRipPerc > 0.05;
temp = sum(temp,1);
i = temp == 0;
suRipPerc(:,i) = nan;
ax.pop.suRipPerc = axes('parent' , F , 'units' , 'inches' , 'position' , [1 axbutt - 0.75 plotWidth 0.5] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8);
hold on;
ax.pop.suRipPerc.XLim = [xaxTickz(1) , xaxTickz(end)];
% ax.pop.suRipFrac.XTick = [0:60:ax.pop.ripEvokedHz.XLim(2)];
ax.pop.suRipPerc.YLim = [0 5];
ax.pop.suRipPerc.YLabel.String = {'% Spikes';'in Ripple'};
for s = 1 : size(suRipPerc , 2)
    plot(ax.pop.suRipPerc , 2.5 + xbins(1:end-1)/60 , suRipPerc(: , s) , 'color' , rcols(s,:) , 'linestyle' , '-' , 'linewidth' , 0.5 , 'marker' , 'none');
end

axbutt = ax.pop.suRipPerc.Position(2);
%% Add another patch of states below for clarity

ax.pop.states2 = axes('parent' , F , 'units' , 'inches' , 'position' , [1 axbutt - 0.2 plotWidth 0.175] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'none' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 8);
hold on;
ax.pop.states2.XLim = [xaxTickz(1) , xaxTickz(end)];
% ax.pop.states2.XTick = [0 : 60 : ripples.expInfo.stopHr*60];
ax.pop.states2.YLim = [0 1];
tic
% Plot states per each second (True line width in inches is 'linewidth'/72
sval = currentdata.state.val;
ssec = currentdata.state.sec;
for s = 1 : length(sval)
    if ~isnan(sval(s))
        t1 = ssec(s);
        patch(ax.pop.states2 , [t1-1,t1,t1,t1-1]/60 , [0 0 1 1] , stateCols(sval(s),:) , ...
            'edgecolor' , 'none' , 'facealpha' , 0.5);
    end
end
toc

axbutt = ax.pop.states2.Position(2);
%% Add x axis at the bottom

ax.pop.xax = axes('parent' , F , 'units' , 'inches' , 'position' , [1 axbutt - 0.1 plotWidth 0.175] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'none' , 'XColor' , 'k' , 'fontname' , 'arial' , 'fontsize' , 8 , 'TickDir'  , 'out');
hold on;
ax.pop.xax.XLim = [xaxTickz(1), xaxTickz(end)];
ax.pop.xax.XTick = xaxTickz;
ax.pop.xax.XLabel.String = 'Minutes';

%% Analyses Plots

% "Ripple-Triggered Averages" with different event types

% (1) Ripple bursts, where the cell fires at least 2 spikes in at least 1
% of the ripples

% Find Ripple Bursts and the start and end of them
% p = currentdata.ripples.peakSec(1:end-1);
% iri = 1000 * (currentdata.ripples.peakSec(2:end) - p);
% 
% a = find(iri<500);
% lI = currentdata.ripples.rangeSec;
% % Find start and end times of each ripple bursts
% for aa = 1 : length(a)
%     ripBurstLimSec(aa,:) = [lI(a(aa),1) , lI(a(aa)+1,2)];
% end

% Only keep the bursts that are totally in NREM

%% Make Figure visible
set(F,'visible' , 'on');
%% single unit rastor plot (v slow maybe just do ripples)
% TL determined this is not such a useful plot..

% ax.suRipRastor = axes('parent' , control_panel , 'units' , 'inches' , 'position' , [1 ax.ripEvokedHz.Position(2) - 5.0 4 1] , 'box' , 'off' , ...
%     'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'k' , 'fontname' , 'arial' , 'fontsize' , 8 , 'XLim' , [0 6*60] , 'XTick' , [0:60:6*60] , 'TickLength' , [0 0]);
% hold on;
% ax.suRipRastor.YLabel.String = 'Unit';
% ax.suRipRastor.YLim = [0 length(spikes.times) - 1];
% ax.suRipRastor.XLim = [0 , ripples.expInfo.stopHr*60];
% ax.suRipRastor.XTick = [0:60:ax.ripEvokedHz.XLim(2)];
%
% lI = currentData.rangeSec;
% for s = 1 : length(spikes.times)
%     spt = spikes.times{s};
%     ripEvoked = [];
%     for i = 1 : size(lI,1)
%         temp = spt(spt >= lI(i,1) & spt < lI(i,2));
%         ripEvoked = [ripEvoked; temp];
%         clear temp;
%     end
%     plot(ax.suRipRastor , repmat(ripEvoked'/60,2,1) , repmat([s-1,s]',1,length(ripEvoked)) , 'linestyle' , '-' , 'color' , rcols(s,:) , 'marker' , 'none' , 'linewidth' ,0.5);
% end

%% Callback functionz and other nested functions

% function cb_states(~ , evt)
% global stateList ripInclude data;
% type = evt.Source.Tag;
% val = evt.Source.Value;
%
% state = find(strcmp(stateList , type));
% ripInclude(data.ripList.states == state) = val;
% clear type val state;

% function cb_rescale(~ , ~)
%
% function cb_sort(~ , evt)
% global h data ripInclude;
% temp = fieldnames(h.sort);
% % set all colors to grey and values to zero
% for t = 1 : length(temp)
%     h.sort.(temp{t}).BackgroundColor = [0.94 0.94 0.94];
%     h.sort.(temp{t}).UserData = 0;
% end
% clear temp t;
% h.sort.(evt.Source.Tag).BackgroundColor = [0.94 0.94 0.5];
% h.sort.(evt.Source.Tag).UserData = 1;
%
% function cb_ripplyzeIt(~ , evt)
% global h data ripInclude stateList control_panel;
%
% for n = 1 : 3
%     h.(evt.Source.Tag).BackgroundColor = [0.5 0.94 0.94];
%     pause(0.1);
%     h.(evt.Source.Tag).BackgroundColor = [0.94 0.94 0.94];
%     pause(0.1);
% end
% % restore logical indices to include all the data
% ripInclude = logical(ones(length(data.ripList.peakPosition) , 1));
%
% % Only include the selected states
% hfn = fieldnames(h.state);
% for hh = 1 : length(hfn)
%     type = h.state.(hfn{hh}).Tag;
%     val = h.state.(hfn{hh}).Value;
%     state = find(strcmp(stateList , type));
%     ripInclude(data.ripList.states == state) = val;
% end
%
% % Now reorder ripInclude and the actual data sets
% % based on value of h.descend and user-selected category
% if h.descend.Value == 1
%     direction = 'descend';
% else
%     direction = 'ascend';
% end
% % Find which sorting variable has been chosen
% % h.sort
% temp = structfun(@(z) z.UserData == 1 , h.sort , 'uniformoutput' , false);
% fn = fieldnames(temp);
% categoryIndx = cell2mat(struct2cell(temp));
% [~ , order] = sort(data.ripList.(fn{categoryIndx})(ripInclude) , direction);
%
% fn = fieldnames(data.ripList);
% for f = 1 : length(data.fn)
%     temp = data.ripList.(fn{f})(ripInclude);
%     currentData.(fn{f}) = temp(order);
% end
% setappdata(control_panel , 'currentData' , currentData);
%
% h.index.String = '1';
% % Clear all axes
% makeRipplePlots(currentData , str2num(h.index.String));
%
% function cb_plus(~ , ~)
% global h control_panel
%
% currentData = getappdata(control_panel , 'currentData');
%
% % Change h.index.Value (or string)
% val = str2num(h.index.String);
% fn = fieldnames(currentData);
% if val ~= size(currentData.(fn{1}),1)
%     h.index.String = num2str(str2num(h.index.String) + 1);
% end
%
% % Run Plotting function
% makeRipplePlots(currentData , str2num(h.index.String));
%
% function cb_minus(~ , ~)
% global h control_panel
%
% currentData = getappdata(control_panel , 'currentData');
%
% % Change h.index.Value (or string)
% val = str2num(h.index.String);
% if val ~= 1
%     h.index.String = num2str(val - 1);
% end
%
% % Run Plotting function
% makeRipplePlots(currentData , str2num(h.index.String));
%
% % function makeRipplePlots(d , indx)
% % % clear all axes..
% %
% % %
% % %    c = get(ax , 'children');
% % %             for cc = 1 : length(c)
% % %                 delete(c(cc));
% % %             end
%
%


function [limz] = findLims(cellArray)
limz = [min(cell2mat(cellfun(@(z) min(z) , cellArray , 'uniformoutput' , false))) , max(cell2mat(cellfun(@(z) max(z) , cellArray , 'uniformoutput' , false)))];

% function [butt] = scaleAx(~ , evt)
% global ax butt
%
% % left click
% if evt.Button == 1
%     butt = evt.IntersectionPoint(1);
% end
%
% % right click
% if evt.Button == 3
%     if exist('butt') == 1
%         if ~isempty(butt)
%             if evt.IntersectionPoint(2) > butt(1)
%                 butt(2) = evt.IntersectionPoint(2);
%                 fn = fieldnames(ax.pop)
%                 for f = 1 : length(fn)
%                     d = round(diff(butt)/2 + butt(1),0);
%                     ax.pop.(fn{f}).XLim = butt;
%                 end
%                 ax.pop.xax.XTick = [butt(1) , d , butt(2)];
%             end
%         end
%     end
%     butt = [];
% end
%
% function  cb_xaxTickz(~,~)
% global ax xaxTickz
% fn = fieldnames(ax.pop)
% for f = 1 : length(fn)
%     ax.pop.(fn{f}).XLim = [xaxTickz(1) , xaxTickz(end)];
% end
% ax.pop.xax.XTick = xaxTickz;





