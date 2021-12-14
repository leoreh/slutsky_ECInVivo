function [ripDetect , ripCount] = TL_ripDetect_tomarMcHugh2020(lfp , sleep , mua)

%% Tomar et al 2020:
% SPW-R events were detected using modifications to the method described in
% (Csicsvari et al., 1999b). Previously selected LFP channels were first
% band-pass filtered (80-250 Hz) using a 69 order Kaiser-window FIR zero-phase shift filter.
% Subsequently, the absolute value of Hilbert transform (instantaneous ripple power)
% was smoothed with 50 ms Gaussian window and candidate SWR-events were detected
% as periods where magnitude exceeded 3 standard deviation (SD) above its
% mean for longer than 30 ms. Further, the initiation and termination points
% of candidate SWR events were defined as points when the magnitude returned
% to the mean. Summed multi-unit activity (MUA) across all pyramidal cells fired
% during any given recorded session was converted to instantaneous firing rate
% using time bin size equal to the underlying LFP trace sampling rate and smoothed,
% to allow detection of firing bursts using the same thresholds as described for
% the candidate SPW-R events detection. Candidate SPW-R events not coincident
% with MUA bursts were excluded from subsequent analysis.

%% do math


% Construct coefficients for butterworth bandpass filter
[c1 c2] = butter(5,[80 250]/(1250/2),'bandpass');

% Run butterworth Filter
f = filtfilt(c1,c2,double(lfp));
clear c1 c2;

% Run Hilbert Transform
h = hilbert(f);

% Extract transformed signal and phase from Hilbert Transform
hilb = h;
amp = abs(h);
phas = angle(h);
% clear zlfp f h;

% Smooth with moving gaussian filter
sm = smoothdata(amp , 'gaussian' , round(50/(1000/1250),0));

logNREM = sleep == 4;
logNREM = repelem(logNREM , 1250);
logNREM = logNREM(1:length(lfp));

mnsm = mean(sm(logNREM));
sdsm = std(sm(logNREM));

zsm = (sm - mnsm)/sdsm;

%% First pass, find ripples based on passing z = 0

thresholded = zsm > 0;

% [start] : index for start of each ripple
% [stop] : index for stop of each ripple
start = find(diff(thresholded)>0);
stop = find(diff(thresholded)<0);
% Exclude last ripple if it is incomplete
if length(stop) == length(start)-1
    start = start(1:end-1);
end
% Exclude first ripple if it is incomplete
if length(stop)-1 == length(start)
    stop = stop(2:end);
end
% Correct special case when both first and last ripples are incomplete
if start(1) > stop(1)
    stop(1) = [];
    start(end) = [];
end

% [ripDetect.lfpIndx] : consolidates start and stop, used moving forward;
ripDetect.lfpIndx = [start,stop];

ripCount.onset = size(ripDetect.lfpIndx,1);

disp(['After z=0 threshold: ' num2str(size(ripDetect.lfpIndx,1)) ' events.']);

%% Second pass, SD threshold
threshOn = 3;
minthreshSamps = 1250*(30/1000);

Li = ripDetect.lfpIndx;
tempStore = [];
for L = 1 : size(Li , 1)
    thresholded = zsm(Li(L,1):Li(L,2)) > threshOn;
    f = find(diff([0;thresholded;0]==1));
    p = f(1:2:end-1);  % Start indices
    y = f(2:2:end)-p;  % Consecutive ones’ counts
    if max(y) >= minthreshSamps
        tempStore(end+1,:) = Li(L,:);
    end
    clear thresholded f p y;
end

ripDetect.lfpIndx = tempStore;
clear tempStore;

ripCount.thresh = size(ripDetect.lfpIndx,1);

disp(['After peak threshold: ' num2str(size(ripDetect.lfpIndx,1)) ' events.']);

%% Throw out ripples that are too short

% Discard ripples that are too long or too short
dur = 1000*(ripDetect.lfpIndx(:,2) - ripDetect.lfpIndx(:,1))/1250;

incl = dur > 30;%

ripDetect.lfpIndx = ripDetect.lfpIndx(incl,:);
ripCount.duration = sum(incl);

% Discard ripples that occur too soon to include a baseline period or too
% late to include a follow up period
incl = ripDetect.lfpIndx(:,1) > 1250*(30/1000) & ripDetect.lfpIndx(:,2) < length(lfp) - 1250*(30/1000);

ripDetect.lfpIndx = ripDetect.lfpIndx(incl,:);
ripCount.edges = sum(incl);

%% Third pass, multi-unit bursting

if ~isempty(mua)
% Bin mua spikes and convert to hz
e = [0 : 10/1000 : length(lfp)/1250];
if e(end) ~= length(lfp)/1250
    e(end+1) = length(lfp/1250);
end
if mua(end) > length(lfp)
    e(end+1) = mua(end);
end

muaBinnedCounts = histcounts(mua , e);
muaBinnedHz = muaBinnedCounts / (10/1000);
if e(end) == mua(end)
    muaBinnedHz(end) = muaBinnedCounts(end)/(e(end)-e(end-1));
end

smMuaBinnedHz = smoothdata(muaBinnedHz , 'gaussian' , 5);
tstamp = e + 5/1000;

m = mean(smMuaBinnedHz);
sd = std(smMuaBinnedHz);
zMua = (smMuaBinnedHz - m)/sd;

% now loop thru ripples and make sure instantaneous firing rate went above
% 3 sds during the ripple

Lt = ripDetect.lfpIndx / 1250;

incl = nan(size(Lt , 1),1);
for L = 1 : size(Lt , 1)
    incl(L) = sum(zMua(tstamp > Lt(L,1) & tstamp < Lt(L,2)) > threshOn) > 0;
end

ripDetect.lfpIndx = ripDetect.lfpIndx(logical(incl) , :);
ripCount.mua = sum(incl);
end



