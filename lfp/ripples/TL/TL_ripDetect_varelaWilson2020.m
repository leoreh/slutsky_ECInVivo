function [ripDetect , ripCount] = TL_ripDetect_varelaWilson2020(lfp , sleep , mua)

%% Info

% Ripple detection BASED on Valera & Wilson 2020

%% Varela & Wilson 2020
% The SWR detection algorithm detected times when the squared, filtered LFP (100–275 Hz)
% had an amplitude above the mean plus three standard deviations for at least 20 ms (mean 
% and standard deviation calculated for the LFP when the animal was quiet). If two SWRs were closer than 20
% ms they were considered a single ripple event. The SWR timestamp was selected as the time with
% the largest absolute value in the ripple filtered LFP.

% ..."and four (KCs, SWRs) standard deviations above the average LFP voltage across sleep periods in
% each recording session"

%% Do math

logNREM = sleep == 4;
logNREM = repelem(logNREM , 1250);
logNREM = logNREM(1:length(lfp));

% Construct coefficients for butterworth bandpass filter
[c1 c2] = butter(5,[100 275]/(1250/2),'bandpass');

% Run butterworth Filter
f = filtfilt(c1,c2,double(lfp));
clear c1 c2;

% square the signal
f2 = f.^2;

% Z-score the entire LFP signal to the NREM mean and standard deviation
mn = mean(f2(logNREM));
sd = std(f2(logNREM));

f2z = (f2 - mn)/sd;

%% First pass using sd threshold

thresholded = f2z > 3;

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

disp(['After onset threshold: ' num2str(size(ripDetect.lfpIndx,1)) ' events.']);

%% Merge based on inter-ripple interval

minInterRippleSamples = 1250*(20/1000);

Li = ripDetect.lfpIndx;

tempStore = [];
ripple = Li(1,:);

for i = 2 : size(Li,1)
    if Li(i,1) - ripple(2) < minInterRippleSamples
        % Merge
        ripple = [ripple(1) Li(i,2)]; % TL now expand the end of the first ripple (combine the ripples)
    else
        tempStore = [tempStore ; ripple];
        ripple = Li(i,:);
    end
end
tempStore = [tempStore ; ripple];
if isempty(tempStore)
    disp('Ripple merge failed');
    return
else
    disp(['After ripple merge: ' num2str(length(tempStore)) ' events.']);
    ripCount.merging = size(tempStore , 1);
end

ripDetect.lfpIndx = tempStore;
clear tempStore;

%% Throw out ripples that are too short

% Discard ripples that are too long or too short
dur = 1000*(ripDetect.lfpIndx(:,2) - ripDetect.lfpIndx(:,1))/1250;

incl = dur > 30;%

ripDetect.lfpIndx = ripDetect.lfpIndx(incl,:);
ripCount.duration = sum(incl);

% Discard ripples that occur too soon to include a baseline period or too
% late to include a follow up period
incl = ripDetect.lfpIndx(:,1) > 1250*(30/1000) & ripDetect.lfpIndx(:,2) < length(f2z) - 1250*(30/1000);

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
    incl(L) = sum(zMua(tstamp > Lt(L,1) & tstamp < Lt(L,2)) > 3) > 0;
end

ripDetect.lfpIndx = ripDetect.lfpIndx(logical(incl) , :);
ripCount.mua = sum(incl);
end

