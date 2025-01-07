function [ripDetect , ripCount] = TL_ripDetect_booneFoster2018(allLFP , vars , sleep , mua)

%% Info

% Ripple detection BASED on Boone..Foster 2018. does not use '3ms' min
% duration (this is only 3 samples...). Uses interripple interval and then
% uses min duration (30ms).
%% Inputs

% [allLFP] : cell array, each cell is the lfp.data vector for a given
% channel, should include 2-3 of the best ripple channels (calculated from
% 'TL_rippleBandRMS'
% [vars] : matrix 1 x 4, inputs correspond to peak threshold, onset
% treshold, consecutive ms above peak threshold, and merge threshold in ms
% [sleep] :

%% Boone..Foster 2018
% LFP and spike data were analyzed in these different wake and sleep states.
% SWRs (marked in red in Fig. 2A) were detected as follows: The three electrodes,
% per hemisphere, with the largest amplitude ripple activity were included in the LFP analysis.
% The z-scored LFP signal of each electrode during “slow-wave sleep or quiet wakefulness” (SWS/QW)
% periods was denoised for 60-Hz electrical noise and its 180-Hz harmonic using a second-order IIR
% notch filter. Denoised LFP during SWS/QW periods was filtered in ripple frequency range
% (100–250?Hz) with a fifth-order Butterworth band-pass filter.
% The envelopes of each band-passed LFP trace were the absolute value of its Hilbert transform.
% These envelopes were averaged over the three electrodes/hemisphere and smoothed with Gaussian
% 5-ms standard deviation (SD) smoother. Each event for which the envelope amplitude exceeds
% a threshold of 3 SD above the mean amplitude for more than 3?ms was considered a SWR,
% and SWRs that were less than 20?ms apart were merged and were considered as one extended event.
% The start and end of each SWR were the times when the smoothed envelope crossed its mean value.
% An analysis, in which the parameters of threshold for SWR detection varied from 1 to 6 SD,
% was conducted to determine the robustness of the findings.

% TL, in my hands this estimates more ripples but the ripple duration
% appears normal. Ripple frequency matches Boone..Foster 2018 if I include
% a peak threshold as well.

%% Do math

logNREM = sleep == 4;
logNREM = repelem(logNREM , 1250);
logNREM = logNREM(1:length(allLFP{1}));

summedSmoov = 0;
for a = 1 : length(allLFP)
    lfp = allLFP{a};
    % Z-score the entire LFP signal to the NREM mean and standard deviation
    mn = mean(lfp(logNREM));
    sd = std(lfp(logNREM));
    
    zlfp = (lfp - mn)/sd;
    
    % Construct coefficients for butterworth bandpass filter
    [c1 c2] = butter(5,[100 250]/(1250/2),'bandpass');
    
    % Run butterworth Filter
    f = filtfilt(c1,c2,double(zlfp));
    clear c1 c2;
    
    % Run Hilbert Transform
    h = hilbert(f);
    
    % Extract transformed signal and phase from Hilbert Transform
    hilb{a} = h;
    amp{a} = abs(h);
    phase{a} = angle(h);
    % clear zlfp f h;
    
    % Smooth with moving gaussian filter
    sm{a} = smoothdata(amp{a} , 'gaussian' , round(5/(1000/1250),0));
    summedSmoov = summedSmoov + sm{a};
end

meanSmoov = summedSmoov / length(allLFP);
zSmoov = (meanSmoov - mean(meanSmoov))/std(meanSmoov);

%  TL their protocol says to use 3*sd + mean (eg zscore of 3) but the data were already z-scored,
% and the mean is low anyhow so seems strange to zscore again or to add the
% mean

% threshPeak = 4 * sdRaw + mean(meanSmoov);
threshPeak = vars(1);
threshOn = vars(2);
peakMs = vars(3);
mergeMs = vars(4);

peakSamps = (peakMs/1000)/(1/1250);
mergeSamps = (mergeMs/1000)/(1/1250);
%% First pass, find ripples based on onset threshold

thresholded = zSmoov > threshPeak;

pos = find(diff(thresholded) == 1) + 1;
neg = find(diff(thresholded) == -1);
if length(pos) > length(neg)
    pos = pos(1:end-1);
else if length(neg) > length(pos)
        neg = neg(2:end);
    end
end
tempStore = [pos , neg];
tempStore = tempStore(neg-pos>=peakSamps , :);

ripDetect.lfpIndx = tempStore;
ripCount.onset = size(tempStore , 1);
clear tempStore thresholded pos neg;
%% Merge based on min interripple interval

% Li = ripDetect.lfpIndx;
% tempStore = [];
% currRip = ripDetect.lfpIndx(1 , :);
% for i = 2 : size(Li , 1)
%     if Li(i,1) - currRip(2) < mergeSamps
%         currRip = [currRip(1) , Li(i,2)];
%     else 
%         tempStore = [tempStore ; currRip];
%         currRip = Li(i , :);
%     end
% end
% 
% ripDetect.lfpIndx = tempStore;
% ripCount.IRI = size(tempStore , 1);
% clear tmepStore currRip;

%% Expand ripple boundaries to z(abs(hilb)) == threshOn
thresholded = zSmoov > threshOn;

pos = find(diff(thresholded) == 1) + 1;
neg = find(diff(thresholded) == -1);

Li = ripDetect.lfpIndx;
tempStore = [];
currRip = Li(1 , :);
for L = 2 : size(Li)
    
    p = pos(pos < currRip(1));
    p = p(end);
    n = neg(neg > currRip(2));
    n = n(1);
    
    if ~isempty(p) & ~isempty(n)
        currRip = [p , n];
    end
    if Li(L,1) <= currRip(2)
        currRip = [p , Li(L,2)];
    else
        tempStore = [tempStore; currRip];
        currRip = Li(L , :);
    end
    
end

ripDetect.lfpIndx = tempStore;
ripCount.onset = size(tempStore , 1);
clear tempStore p n pos neg Li L;

% %% Throw out ripples that are too long or too short or have artifact in raw tdt trace
% 
% % Discard ripples that are too long or too short
% dur = 1000*(ripDetect.lfpIndx(:,2) - ripDetect.lfpIndx(:,1))/1250;
% 
% incl = dur > 30;%
% 
% ripDetect.lfpIndx = ripDetect.lfpIndx(incl,:);
% ripCount.duration = sum(incl);
% 
% % Discard ripples that occur too soon to include a baseline period or too
% % late to include a follow up period
% incl = ripDetect.lfpIndx(:,1) > 1250*(30/1000) & ripDetect.lfpIndx(:,2) < length(meanSmoov) - 1250*(30/1000);
% 
% ripDetect.lfpIndx = ripDetect.lfpIndx(incl,:);
% ripCount.edges = sum(incl);

%% Use MUA criteria if user decides

if ~isempty(mua)
    lfp = allLFP{1};
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
        incl(L) = sum(zMua(tstamp > Lt(L,1) & tstamp < Lt(L,2)) >= 3) > 0;
    end
    
    ripDetect.lfpIndx = ripDetect.lfpIndx(logical(incl) , :);
    ripCount.mua = sum(incl);
end



