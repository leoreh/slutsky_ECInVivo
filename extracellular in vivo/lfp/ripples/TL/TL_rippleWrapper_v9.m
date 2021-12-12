function [ripples] = TL_rippleWrapper_v9(pathname , set , rms , saveD)

% v7: edited to work with TL new spikewrapper_v2 +
%v6: in case of no sleep data, calculates theta/delta ratio and sorts into
%nrem, wake, and unknown.

%% Inputs:

% [pathname] : string pathname to folder containing recording data
% (eg. 'E:\refaela\wt2\200820');
% [set] : set number to analyze
% [rms] : 1 to use rms calculation for ripple channels, 0 to just pick the
% 1st electrode on each tetrode, ranked by most multiunit spiking on 1st
% baseline day (simplest..)
% [saveD] : 1 to save the ripples structure, 0 to not save
%% Ripple Detection Parameters from Literature

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

% fix pathname string
if ~strcmp(pathname(end) , filesep)
    pathname = [pathname , filesep];
end

% this is not really necessary can adjust later..
% Load expInfo from original date
temp = dir([pathname, '*Raw' , '*expInfo.mat']);
% i = find(cell2mat(cellfun(@(z) ~isempty(z) , strfind({temp.name} , 'Raw') , 'uniformoutput' , false)));
load([pathname , temp.name]);

setI = find(cell2mat(expInfo.set) == set , 1);
dt = expInfo.date{setI};
fs = strfind(pathname , filesep);
saveFol = [pathname(1:fs(end-1)) , num2str(dt) , filesep];

%% Load expInfo from the excel spreadsheet containing experiment info..

% Check directory in pathname ---------------------------------------------
% [fname] : names of files in pathname, including their suffix (eg .mat ,.xls)
temp = dir(pathname);
fnames = {temp(~[temp.isdir]).name}';
clear temp;

% Load excel spreadsheet with experiment and analysis parameters and
% extract info ------------------------------------------------------------

% Load expInfo
temp = dir([saveFol, '*Raw' , '*expInfo.mat']);
% i = find(cell2mat(cellfun(@(z) ~isempty(z) , strfind({temp.name} , 'Raw') , 'uniformoutput' , false)));
load([saveFol , temp.name]);
% load([pathname , temp(i).name]);
clear temp i;

setI = find(cell2mat(expInfo.set) == set,1);

% load emgInfo
d = dir([expInfo.saveFol{setI} , filesep , '*EMG' , '*expInfo.mat']);
emgInfo = load([expInfo.saveFol{setI} , filesep , d.name]);
emgInfo = emgInfo.expInfo;
clear d;

% directory to emg.dat
 d = dir([expInfo.saveFol{setI} , filesep , '*EMG' , '*.dat']);
 emgDir = [expInfo.saveFol{setI} , filesep , d.name];
 clear d;
%% First check if 'best' ripple channels have been detected yet by calculating RMS
% If not, run TL_rippleBandRMS
% If so, load 'RMS.mat' from subject's path (eg. 'E:\refaela\wt2\')

f = strfind(pathname , filesep);
f = f(end-1);

subjPath = pathname(1:f); clear f;

if rms
    d = dir([subjPath , '*beztRipChans*']);
    if ~isempty(d)
        load([subjPath , filesep , 'beztRipChans.mat']);
    else
        [beztRipChans] = TL_rippleBandRMS_v2(subjPath , 1);
        %     [beztRipChans] = TL_rippleBandRMS(subjPath , baseline , sleep , saveLFP , saveRMS);
    end
    clear subjPath;
    
    % Now find best channel on each of the 4 tetrodes, and take the 3 best
    % of these
    allCh = [1:16];
    ut = unique(beztRipChans.tet);
    for t = 1 : length(ut)
        currCh = allCh(beztRipChans.tet == ut(t));
        [~ , r] = ismember(beztRipChans.rank , currCh);
        i = find(r>0);
        ripChans(t) = currCh(r(i(1)));
    end
    [~ , r] = ismember(beztRipChans.rank , ripChans);
    i = find(r>0);
    ripChans = ripChans(r(i(1:3)));
    
else
    
    % v dirty way of doing this
    d = dir(subjPath);
    d = {d([d.isdir]).name};
    i = cell2mat(cellfun(@(z) ~strcmp(z(1) , '.') , d , 'uniformoutput' , false));
    recFols = [d(i)];
    r = sort(recFols); r = r{1};
    t = dir([subjPath , r , filesep , '*spktimes.mat']);
    load([t.folder , filesep , t.name]);
    [~ , tets] = sort(cell2mat(cellfun(@(z) length(z) , spktimes , 'uniformoutput' , false)) , 'descend');
    ripChans = 1 + (tets-1)*4;
    ripChans(ripChans>length(expInfo.channels{1})) = length(expInfo.channels{1});
    clear d i recFols r t;
    
end

%% Extract sleep from the recording folder
temp = dir([pathname '*sleep.mat']);
ch2use = ripChans(1);
if isempty(temp)
    sleep = TL_extractSleep(pathname , expInfo , ch2use);
    save([pathname , 'sleep.mat'] , 'sleep');
else
    load([pathname , 'sleep.mat']);
end
% sleep.mat contains the variable 'labels'
%% Load multi unit spike data

temp = fnames{cell2mat(cellfun(@(z) ~isempty(strfind(z,'spktimes')) , fnames , 'uniformoutput' , false))};
muSpkSec = load([pathname , temp]);
% [r,~] = find(expInfo.tetrodes == ripChans(1));

% some mu spike files refaela sent were as integer indices..others as
% time..so check format and convert to seconds if necessary
tUse = tets(1);
if sum(~rem(muSpkSec.spktimes{tUse},1)) == length(muSpkSec.spktimes{tUse})
    muSpkSec = muSpkSec.spktimes{tUse}/expInfo.fs{setI};
else
    muSpkSec = muSpkSec.spktimes{tUse};
end
expInfo.ripChan{1} = ripChans(1);
clear temp r;

%% Load single unit spike data
% temp = fnames{cell2mat(cellfun(@(z) ~isempty(strfind(z,'spikes.cellinfo')) , fnames , 'uniformoutput' , false))};
% temp = load([pathname , temp]);
% suSpkSec = temp.spikes.times;
% clear temp;

%%
% calc total time of recording
% dv = cellfun(@(z) datevec(z)' , expInfo.duration , 'uniformoutput' , false);
% dv = cell2mat(dv);
% dv = sum(dv , 2);
% dv = dv(4:end);
% totHr = dv(1) + dv(2)/60 + dv(3)/360;

% for some reason this calculation of time (and nsec) is coming out differently than
% the output length of getLfp (missing a few hours in 1 case 210611 ...
totHr = sum(cell2mat(expInfo.blockduration))/3600;
segments = [0:12:totHr];
% TL temporary fix
if expInfo.EndSec{setI} > 0
    segments(end) = expInfo.EndSec{setI}/3600;
end
if segments(end) ~= totHr
    segments(end+1) = totHr;
end
% segments = [0,12,22];
offset = 0;

for e = 1 : length(segments) - 1
    
    strtSec = segments(e)*60*60;
    if strcmp(segments(e+1),'inf')
        stopSec = inf;
    else
        stopSec = segments(e+1)*60*60; %% problem here reading 2 digit integer..
    end
    
    for r = 1 : length(ripChans)
        [filtLfp{r} , ~] = TL_getLFPtest('basepath' , saveFol , 'type' , 'filt' , 'nchans' , length(expInfo.channels{setI}) , ...
            'ch' , ripChans(r) , 'fs' , 1250 ,  'extension' , 'lfp' , 'savevar' , false , 'forceL' , true , 'interval' , [strtSec stopSec]);
    end
    
    % chunk data into smaller bits due to higher sampling rate of
    % unfiltered trace
    if stopSec ~= inf
        ufRange = strtSec:60*60:stopSec;
    else
        ufRange = strtSec:60*60:length(labels);
%         ufRange(end) = inf;
    end
    ufRange(end) = inf;
    
    unfiltLfp = [];
    unfiltTs = [];
    for u = 1 : length(ufRange)-1
        % Needs smaller chunk sizes (divide by 10)
        [temp, ~] = TL_getLFPtest('basepath' , pathname , 'type' , 'unfilt' , 'nchans' , length(expInfo.channels{setI}) , ...
            'ch' , ripChans(1) , 'fs' , 12500 , 'extension' , 'lfp' , 'savevar' , false , 'forceL' , true , 'interval' , [ufRange(u) , ufRange(u+1)] , 'cf' , []);
        unfiltLfp = [unfiltLfp ; temp.data];
        unfiltTs = [unfiltTs ; temp.timestamps];
        
        clear temp;
    end
    clear r strtHr stopHr;
      
    recLength(e) = length(filtLfp{1})/1250;
    
    if stopSec > length(labels)
        currSleep = labels([strtSec + 1 : end]);
    else
        currSleep = labels([strtSec + 1 : stopSec]);
    end
    
    



 % Needs smaller chunk sizes (divide by 10)

 
 % TL now automatically just using channel 1
  emgData = double(bz_LoadBinary(emgDir, 'duration', stopSec - strtSec , ...
        'frequency', emgInfo.fs{setI} , 'nchannels', 4, 'start', strtSec,  ...
        'channels', 1));

    [d , c] = TL_ripDetect_v4(filtLfp , unfiltLfp , emgData , [4,1,10,20] , currSleep , muSpkSec);
    d.lfpIndx = d.lfpIndx + offset;
    Lth(e) = size(d.lfpIndx,1);
    tempDetect{e} = d;
    tempCount{e} = c;
    offset = offset + length(filtLfp{1}.data);
    clear d c lfp filtLfp unfiltLfp;
end

expInfo.recLengthSec = sum(recLength);
ripples.expInfo = expInfo;

ffn = fieldnames(tempDetect{1});
for f = 1 : length(ffn)
    tempD = [];
    for e = 1 : length(tempDetect)
        tempD = [tempD; tempDetect{e}.(ffn{f})];
    end
    ripples.detect.(ffn{f}) = tempD;
    clear tempD;
end

clear f m;

ffn = fieldnames(tempCount{1});

for f = 1 : length(ffn)
    tempC = [];
    for e = 1 : length(tempCount)
        tempC = tempC + tempCount{e}.(ffn{f});
    end
    ripples.count.(ffn{f}) = tempC;
    clear tempC
end

ripples.set = set;

clear f m;

clear fn;

if saveD
    save([pathname , 'ripples_' num2str(set) '.mat'] , 'ripples')
end
