function [ripDetect , ripCount] = TL_ripDetect_v4(filtLFP , unfiltLFP , emgData , vars , sleep , mua)

%% Ripple detection based on Boone..Foster 2018 with changes to onset threshold z-scored magnitude
% minimum inter-ripple interval, and min time above peak z-scored
% threshold. Parameters were tested using TL_compareRipDetectMethods for
% optimization. Does not use multi-unit activity criteria but this can be
% saved and added to the ripple structure separately. This code just gives
% the onset and offset indices of each ripple, and displays how many
% ripples were lost at each filtering step

% v2 : just clearing some unneeded variables to avoid storage problems
% v3 : chunking the data during hilbert to bins of 12 hours to avoid storage problems..


% INPUTS:
% [filtLFP] : cell array, each cell is the LFP low pass filtered at 450 and
% sampled at 1250 hz at the few channels with the best ripples (recommend
% 3 ch's). Previously outputted fro TL_getLFPtest
% [unfiltLFP] : unfilteredLFP sampled at 12500 and no filtering, which will
% be used to visual the spiking during the ripples and the raw ripple
% frequency
% [emgData] : vector of emg data same time span as filtLFP. Frequency will
% be back calculated later
% [vars] : 4 value vector of variables used for ripple detection: 
% [(1) z-scored hilb transform amplitude used as peak threshold detection
% criterion, (2) z-scored "" for onset threshold, (3) ms that that peak
% threshold must surpassed, (4) ms threshold to merge temporally adjacent
% ripples] .... [4,1,10,20] works well, can play around from these values

%% Inputs:
% [allLFP] : cell array, each cell is the lfp.data output from TL_getLFP
% for the 3 ripple channels chosen (find best rms channel from each
% tetrode, then take the 3 best of these)

%% Original Boone..Foster 2018 detection protocol this is based on

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

%% Other Ripple Detection Parameters from Literature

% Tomar et al 2020:
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

% Roux..Buzsaki et al 2017
% The root-mean-square (RMS) of the bandpass-filtered (80–250-Hz) signal
% was computed in two running windows, long (2 s; RMS1) and short (8 ms; RMS2).
% Ripples were defined as events with RMS2 exceeding 3ª RMS1 (range: 3–3.5) for at least 8 ms
% (ref. 31).

% Norimoto et al 2018
% The root mean square (RMS) of the bandpass filtered (140–250 Hz) signal was
% computed in two running windows, long (2 s; RMS1) and short (8 ms; RMS2). Ripples
% were defined as events with an RMS2 exceeding 4ªRMS1 for at least 8 ms

% Varela & Wilson 2020
% The SWR detection algorithm detected times when the squared, filtered LFP (100–275 Hz)
% had an amplitude above the mean plus three standard deviations for at least 20 ms (mean
% and standard deviation calculated for the LFP when the animal was quiet). If two SWRs were closer than 20
% ms they were considered a single ripple event. The SWR timestamp was selected as the time with
% the largest absolute value in the ripple filtered LFP.

% Boone..Foster 2018
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

% Altimus..Foster 2015
% For LFP analysis, 1 representative electrode was selected from each tetrode,
% and the 3 tetrodes with the strongest ripple activity were selected for further analysis.
% The z-score of each channel was taken to normalize for differences in recordings between animals.
% On these 3 tetrodes, the z-score of the LFP was band-pass filtered between 100 and 300 Hz,
% and the absolute value of the Hilbert transform of this filtered signal was then smoothed
% (Gaussian kernel, SD = 10 ms). This processed signal was then averaged across all tetrodes,
% and ripple events were identified as local peaks with an amplitude >5 SD above the mean
% during periods that the mouse's velocity was <1 cm/s. The start and end boundaries of each
% ripple event were defined as the point when the signal crossed mean + 2 SD. Finally,
% to prevent inappropriate splitting of ripple events, those that occurred within 100 ms
% of the end of the previous event were combined into a single event. The same ripple finding
% algorithm was also used for high gamma (frequency range 65-140). For theta power analysis,
% the raw LFP trace was band-pass filtered between 4 and 12 Hz during bouts of activity when
% the mouse's velocity exceeded 4 cm/s. The absolute value of the Hilbert transform for the
% filtered signal was z-scored. The power spectral density was then calculated using Welch's
% power spectral density estimate of this data. Low gamma (frequency range 30-90) was
% examined with the same algorithm used in the theta analysis as it is predominantly expressed during bouts of activity.

% Suh..Tonegawa 2013
% One electrode from each tetrode that had at least one cluster was considered for EEG analysis.
% EEG signal of each electrode was denoised for 60 Hz electric noise and its 180 Hz harmonic
% using a second-order IIR notch filter. Denoised EEG was filtered at ripple frequency range
% (100–240 Hz) with a fifth-order Butterworth band-pass filter. The envelopes of each band-passed
% EEG were obtained using the absolute value of its Hilbert transform and these envelopes were
% averaged over all electrodes. After applying a Gaussian smoother with 5 ms standard deviation,
% the averaged envelope was z-scored. Events that passed 5 standard deviations
% (i.e., mean + 5 SD of averaged non-z-scored envelope) for more than 3 ms were considered as ripples,
% and ripples that were less than 20 ms apart were merged and were considered as one extended ripple.
% The beginning and end of each ripple were considered as where the smoothed envelope
% crossed its mean value (i.e., zero for z-scored signal). Ripples events that happened
% when mice were not immobilized were excluded. Mice were considered as immobilized when
% their head speed was below 0.5 cm/s. Ripple power was obtained by applying Welch’s method
% on each individual non-z-scored nonenveloped ripple and then averaging over calculated powers.

% Cheng & Frank 2008
% On a given day, the tetrode with the largest number of isolated neurons was chosen for ripple analysis.
% The LFP envelope was determined by Hilbert transforming the LFP signal from
% that tetrode after band pass filtering between 150 and 250 Hz. Events were detected
% if the envelope exceeded threshold for at least 15 ms. Events included times around
% the triggering event during which the envelope exceeded the mean.
% Overlapping events were combined. A high threshold of mean + 6 SD of the envelope
% was used for detecting events that would meet most common definitions of ripples.
% HFEs were detected with a lower threshold of mean + 3 SD.

%% Do math

stateList = {'activeWake', 'quietWake' , 'lightSleep' , 'NREM' , 'REM' , 'N/REM' , 'bin' , 'unknown'};

% -------------------------------------------------------------------------
% Find NREM logical indices in filtLFP trace and trim data to equal size of
% sleep vector.

% [logNREM] : index in 1250khz LFP for NREM
nremI = find(strcmp(stateList , 'NREM'));
logNREM = sleep == nremI; clear nremI;
logNREM = repelem(logNREM , 1250); % assumes lfp sampling rate of 1250 hz

% trim data as necessary
L = length(filtLFP{1}.data);
if L < length(logNREM)
    logNREM = logNREM(1:length(filtLFP{1}.data));
end

if L > length(logNREM)
    for a = 1 : length(filtLFP)
        filtLFP{a}.data = filtLFP{a}.data(1:length(logNREM));
        filtLFP{a}.timestamps = filtLFP{a}.timestamps(1:length(logNREM));
%         filtLFP = cellfun(@(z) z(1:length(logNREM)) , filtLFP , 'uniformoutput' , false);
    end
end
clear L;

% -------------------------------------------------------------------------

summedSmoov = 0;
summedPhase = 0;
summedFilt = 0;

chunkSamps = 12*60*60*1250;

% loop thru ripple channels
for a = 1 : length(filtLFP)
    lfp = filtLFP{a}.data;
    
    % Z-score the entire LFP signal to the NREM mean and standard deviation
    mn = mean(lfp(logNREM));
    sd = std(lfp(logNREM));
    zlfp = (lfp - mn)/sd;
    clear lfp; 
    
    % Construct coefficients for butterworth bandpass filter
    [c1 , c2] = butter(5,[100 250]/(1250/2),'bandpass');
    
    % Run butterworth Filter
    f = filtfilt(c1,c2,double(zlfp));
    clear c1 c2 zlfp;
    
    % chunk this section to save space --------------------------
    %     if length(f)
    chunks = 0 : chunkSamps : length(f);
    if chunks(end) ~= length(f)
        if chunks == 0
            chunks = [0,length(f)];
        else
            if length(f) - chunks(end-1) > 15*60*60
                chunks(end+1) = length(f);
            else chunks(end) = length(f);
            end
        end
    end
    h = [];
    sm = [];
    ph = [];
    
    for ch = 1 : length(chunks) - 1
        
        %         if ch == 1
        %             d = f(chunks(ch) : chunks(ch+1));
        %         else
        d = f(chunks(ch) + 1 : chunks(ch+1));
        %         end
        
        % Run Hilbert Transform
        h = hilbert(d);
        %         h = [h; hilbert(d)];
        clear d;
        
        % Extract transformed signal and phase from Hilbert Transform
        %     hilb{a} = h;
        amp = abs(h);
        ph = [ph ; angle(h)];
        clear h;
        
        % Smooth with moving gaussian filter
        sm = [sm ; smoothdata(amp , 'gaussian' , round(5/(1000/1250),0))];
        clear amp;
    end
    summedSmoov = summedSmoov + sm;
    clear sm;
    
    summedPhase = summedPhase + ph;
    
    summedFilt = summedFilt + f;
    clear f;
end

meanSmoov = summedSmoov / length(filtLFP);
clear summedSmoov;

zSmoov = (meanSmoov - mean(meanSmoov))/std(meanSmoov);

meanPhase = summedPhase / length(filtLFP);
clear summedPhase;

meanFilt = summedFilt / length(filtLFP);

zFilt = (meanFilt - mean(meanFilt))/std(meanFilt);

clear summedFilt;

% threshPeak = 4 * sdRaw + mean(meanSmoov);
threshPeak = vars(1);
threshOn = vars(2);
peakMs = vars(3);
mergeMs = vars(4);

peakSamps = (peakMs/1000)/(1/1250);
mergeSamps = (mergeMs/1000)/(1/1250);
%% First pass, find ripples based on peak threshold

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

%% Expand ripple boundaries to z(abs(hilb)) == threshOn
thresholded = zSmoov > threshOn;

pos = find(diff(thresholded) == 1) + 1;
neg = find(diff(thresholded) == -1);

Li = ripDetect.lfpIndx;
tempStore = [];
currRip = Li(1 , :);
for L = 2 : size(Li)
    if L == 1201
    end
    p = pos(pos <= currRip(1));
    p = p(end);
    n = neg(neg >= currRip(2));
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

%% Remove ripples that are too short or too long
% TL finds this helps remove artifacts

d = ripDetect.lfpIndx(:,2) - ripDetect.lfpIndx(:,1);
d = 1000*d/1250;
excl = ~ [d < 20 | d > 200];

ripDetect.lfpIndx = ripDetect.lfpIndx(excl , :);
ripCount.length = size(ripDetect.lfpIndx , 1);
clear excl d;

%% Extract phase, envelope, unfiltered trace, and filtered trace in each ripple

Li = ripDetect.lfpIndx;

bufferLfpSamps = (20/1000)*1250;
bufferLfpMs = 1000 * bufferLfpSamps/1250;
bufferTdtSamps = round(24414.06 * bufferLfpMs/1000 , 0);
bufferTdtMs = 1000 * bufferTdtSamps / 24414.06;
filtTs = filtLFP{1}.timestamps;

for L = 1 : size(Li , 1)
    range(1) = Li(L,1) - bufferLfpSamps;
    range(2) = Li(L,2) + bufferLfpSamps;
    range = range(1) : range(2);
    % In case range goes out of recording limits
    buffPre = bufferLfpSamps - sum(range <= 0);
    buffPost = bufferLfpSamps - sum(range > length(meanSmoov));
    
    
    range = range(range > 0 & range <= length(meanSmoov));
%     range<0
    phase{L,1} = meanPhase(range);
    
%     envelope{L,1} = meanSmoov(range);
      envelope{L,1} = zSmoov(range);
      
    % match timestamps in filtered Range to timestamps in unfiltered range
    % need to make estimation first because searching min through entire
    % unfilTs takes too long. estimation is filtTs*
%     [~ , ufr(1)] = min(abs(unfiltTs - filtTs(range(1))));
%     [~ , ufr(2)] = min(abs(unfiltTs - filtTs(range(end))));
ufr = range(1)*10:range(end)*10;
    unfiltered{L,1} = unfiltLFP(ufr) - mean(unfiltLFP(ufr(1) - bufferTdtSamps : ufr(1)));
    
%     filtered{L,1} = meanFilt(range) - mean(meanFilt(range(1:buffPre)));
    filtered{L,1} = zFilt(range);% mean(meanFilt(range(1:buffPre)));
    
    temp = ones(length(range) , 1);
    temp([1:buffPre, end-buffPost+1:end]) = 0;
    ripLog{L,1} = temp; 
    clear temp ufr buffPre buffPost range;
end

ripDetect.phase = phase;
ripDetect.envelope = envelope;
ripDetect.unfiltered = unfiltered;
ripDetect.filtered = filtered;
ripDetect.ripLog = ripLog;

%% Remove artifact traces

% based on detecting rapid slope
artSize = []
for d = 1 : length(ripDetect.unfiltered)
    u = ripDetect.unfiltered{d};
    L = length(u) - 10;
    for i = 1 : L
        slp(i) = u(i+9) - u(i);
    end
    artSize(d) = min(slp);
    clear slp u L;
end
removeArt = artSize < -750;

fn = fieldnames(ripDetect);
for f = 1 : length(fn)
    ripDetect.(fn{f}) = ripDetect.(fn{f})(~removeArt , :);
end
clear artSize removeArt;
%% Extract z-scored multi-unit activity in each ripple and mu spike times

lfp = filtLFP{1}.data;
binRate = 5/1000;
% Bin mua spikes and convert to hz
e = [0 : binRate : length(lfp)/1250];
if e(end) ~= length(lfp)/1250
    e(end+1) = length(lfp)/1250;
end
if mua(end) > length(lfp)
    e(end+1) = mua(end);
end

muaBinnedCounts = histcounts(mua , e);
muaBinnedHz = muaBinnedCounts / binRate;
if e(end) == mua(end)
    muaBinnedHz(end) = muaBinnedCounts(end)/(e(end)-e(end-1));
end
clear muaBinnedCounts;

smMuaBinnedHz = smoothdata(muaBinnedHz , 'gaussian' , 5);
tstamp = e + binRate/2;
tstamp = tstamp(1:end-1);

m = mean(smMuaBinnedHz);
sd = std(smMuaBinnedHz);
zMua = (smMuaBinnedHz - m)/sd;
clear smMuaBinnedHiz

% for ease of plotting and analysis later on, resample zMua to 1/1250 (lfp
% sample rate)
zMua = resample(zMua , (1/1250)^-1 , binRate^-1);

% now loop thru ripples and make sure instantaneous firing rate went above
% 3 sds during the ripple

Li = ripDetect.lfpIndx;
Lt = Li / 1250;

for L = 1 : size(Lt , 1)
    
     range(1) = Li(L,1) - bufferLfpSamps;
    range(2) = Li(L,2) + bufferLfpSamps;
    range = range(1) : range(2);
    
    % In case range goes out of recording limits
    buffPre = bufferLfpSamps - sum(range <= 0);
    buffPost = bufferLfpSamps - sum(range > length(zMua)); 
     range = range(range > 0 & range <= length(zMua));
     
     zm{L,1} = zMua(range);
     
      temp = ones(length(range) , 1);
    temp([1:buffPre, end-buffPost+1:end]) = 0;
    zmuaLog{L,1} = temp; clear temp;
     
    temp = ones(length(zm{L,1}),1);
    temp([1:bufferLfpSamps, end - bufferLfpSamps + 1 : end]) = 0;
    logZmua{L,1} = temp; clear temp;
    
    i = mua >= Lt(L,1) - bufferTdtMs/1000 & mua < Lt(L,2) + bufferTdtMs/1000;
    spkT{L,1} = mua(i) - Lt(L,1);
    clear i range buffPre buffPost;
end

ripDetect.zMua = zm;
ripDetect.logZmua = logZmua;
ripDetect.spkT = spkT;
clear zm logZmua spkT;
%% Remove ripples that do not evoke multi-unit spiking

Li = ripDetect.lfpIndx;
Lt = Li / 1250;
Lt = Lt(:,2) - Lt(:,1);

for L = 1 : size(Lt , 1)
    sp = ripDetect.spkT{L};
    bl(L) = sum(sp < 0) / (bufferLfpMs/1000);
    evk(L) = sum(sp >=0 & sp < Lt(L)) / Lt(L);
end
excl = evk == 0; % not use evk over bl, just evokes at least 1 spk

fn = fieldnames(ripDetect);
for f = 1 : length(fn)
    ripDetect.(fn{f}) = ripDetect.(fn{f})(~excl , :);
end
clear excl Li Lt sp bl evk;
%% Remove ripples with emg artifact

if ~isempty(emgData)
    
    % params for downsample the EMG 
    EmgFs = 1250*length(emgData)/length(filtLFP{1}.data);
    ds = floor(EmgFs / 1250);
    newFs = EmgFs / ds;
    oldLength = length(emgData);
    newLength = oldLength / ds;
    
        % chunk this section to save space --------------------------
    %     if length(f)
    chunkSamps = 10000000;
    chunks = 0 : chunkSamps : oldLength;
    if chunks(end) ~= oldLength
        if chunks == 0
            chunks = [0,oldLength];
        else
            if length(emgData) - chunks(end-1) > 15*60*60
                chunks(end+1) = oldLength;
            else chunks(end) = oldLength;
            end
        end
    end
    
    dsEmgData = [];
    for c = 1 : length(chunks) - 1
    dsEmgData = [dsEmgData ; downsample(emgData(chunks(c)+1: chunks(c+1)) , ds)];
    end
    zEmg = (dsEmgData - mean(dsEmgData))/std(dsEmgData);
    
    % Construct coefficients for butterworth bandpass filter
    [c1 , c2] = butter(5,[100 250]/(1250/2),'bandpass');
    
    % Run butterworth Filter
    f = filtfilt(c1,c2,double(zEmg));
    clear c1 c2 zlfp;

    h = [];
    sm = [];
    ph = [];
    
    % Rechunk
    
      
        % chunk this section to save space --------------------------
    %     if length(f)
    chunkSamps = 10000000;
    chunks = 0 : chunkSamps : newLength;
    if chunks(end) ~= newLength
        if chunks == 0
            chunks = [0,newLength];
        else
            if length(emgData) - chunks(end-1) > 15*60*60
                chunks(end+1) = newLength;
            else chunks(end) = newLength;
            end
        end
    end
    
    for ch = 1 : length(chunks) - 1
        
        %         if ch == 1
        %             d = f(chunks(ch) : chunks(ch+1));
        %         else
        d = f(chunks(ch) + 1 : chunks(ch+1));
        %         end
        
        % Run Hilbert Transform
        h = hilbert(d);
        %         h = [h; hilbert(d)];
        clear d;
        
        % Extract transformed signal and phase from Hilbert Transform
        %     hilb{a} = h;
        amp = abs(h);
        clear h;
        
        % Smooth with moving gaussian filter
        sm = [sm ; smoothdata(amp , 'gaussian' , round(5/(1000/newFs),0))];
        clear amp;
    end
    
    zsm = (sm-mean(sm))/std(sm);
    
    thresholded = zsm > 6;
    fthresh = find(thresholded);
    
    % convert to lfpIndx
    fthresh * 1250 / newFs;
    
    % loop thru ripples
    Li = ripDetect.lfpIndx;
    flag = zeros(size(Li , 1) , 1);
    for L = 1 : size(Li , 1)
        if sum(fthresh >= Li(L,1) & fthresh <= Li(L,2)) > 0
            flag(L) = 1;
        end
    end
    excl = logical(flag);
    fn = fieldnames(ripDetect);
    for f = 1 : length(fn)
        ripDetect.(fn{f}) = ripDetect.(fn{f})(~excl , :);
    end
end