function events = findLFPevents(varargin)

% detects events from LFP signal. designed for burst-suppression and
% ripples. based on bz_FindRipples. main differences are (1) moving average
% implemented via movmean, (2) detects peak after elimination, (3) selects
% best channel automatically, (4) uses only single threshold. more
% information is availble in PreProcessing.docx
% 
% INPUT
%   lfp         struct (see getLFP.m)
%   preset      predefined events {ripples} or (bs)
%   passband    to filter LFP. if empty will use raw LFP
%   winLength   for moving mean 
%   minDur      minimum event duration [ms]
%   maxDur      maximum event duration [ms]
%   interDur    minimum time between events [ms]
%   ch          channel to detect events. if empty will select by RMS, 
%               if > 1 will average signal across specified channels
%   noise       channel to remove artifacts {[]}
%   emgThr      threshold for EMG exclusion {0.9}
%   interval    interval used for standadization {[0 Inf]}
%   basepath    recording session path {pwd}
%   graphics    plot figure {1}.
%   saveVar     save variable {1}.
% 
% OUTPUT
%   events         struct with fields:
%       stamps      n x 2 mat where n is the number of events and the
%                   columns represent the start and end of each event [samps]
%       peaks       the timestamp for each peak [samps] 
%       peakPower   the voltage at each peak [uV]
%       noise       events execluded based on noise channel
%       <params>    as in input
% 
% CALLS
%   fastrms
% 
% TO DO LIST:
%   compare with original script on various data
%   track and save exclusions
%   add raw lfp to graphics
% 
% 01 may 19 LH. updates:
% 07 oct 19 LH  adapted for BS
% 10 oct 19 LH  combined BS and ripples
% 16 oct 19 LH  changed output to samples instead of 
% 18 nov 19 LH  adapted for IIS


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'lfp', @isstruct)
addParameter(p, 'preset', '', @ischar)
addParameter(p, 'passband', [], @isnumeric)
addParameter(p, 'winLength', [], @isnumeric)
addParameter(p, 'minDur', [], @isnumeric)
addParameter(p, 'maxDur', [], @isnumeric)
addParameter(p, 'interDur', [], @isnumeric)
addParameter(p, 'ch', [], @isnumeric)
addParameter(p, 'noise', [], @isnumeric)
addParameter(p, 'emgThr', 0.9, @isnumeric);
addParameter(p, 'interval', [0 Inf], @isnumeric)
addParameter(p, 'basepath', pwd, @isstr);
addParameter(p, 'graphics', true, @islogical)
addParameter(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
lfp = p.Results.lfp;
preset = p.Results.preset;
passband = p.Results.passband;
winLength = p.Results.winLength;
minDur = p.Results.minDur;
maxDur = p.Results.maxDur;
interDur = p.Results.interDur;
ch = p.Results.ch;
noise = p.Results.noise;
emgThr = p.Results.emgThr;
interval = p.Results.interval;
basepath = p.Results.basepath;
graphics = p.Results.graphics;
saveVar = p.Results.saveVar;

% check interval
inter = interval * lfp.fs;
if inter(2) == Inf
    inter(2) = length(lfp.data);
end
if inter(1) == 0
    inter(1) = 1;
end
inter = inter(1) : inter(2);

% prepare ch
if isempty(ch)
    ch = 1 : size(lfp.data, 2);
end

% prepare output
events = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
switch preset
    case 'ripples'
        passband = [130 200];           % [Hz]
        order = 4;
        type = 'butter';
        winLength = 11;                 % [samples]
        
        thr = [2 5];                    % low and high thresholding value [z-score]
        minDur = 0.02;                  % minimum event duration [s]
        maxDur = 0.1;                   % maximum event duration [s]
        interDur = 0.03 * lfp.fs;       % minimum time between events [s]
        
    case 'bs'
        passband = [5 250];
        order = 4;
        type = 'butter';
        winLength = 30;
        
        minDur = 0.05;
        maxDur = 4000;
        interDur = 0.2 * lfp.fs;
        thr = [0.01 1];
        
    case 'iis'
        passband = [1 Inf];
        order = 4;
        type = 'cheby2';
        winLength = 30;
        
        minDur = 0.004;
        maxDur = 4000;
        interDur = 0.2 * lfp.fs;
        thr = [2 1];
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filter    
if ~isempty(passband)
    sig = filterLFP(double(lfp.data(:, ch)), 'type', type, 'fs', lfp.fs,...
        'passband', passband, 'order', order, 'graphics', false);
else
    sig = double(lfp.data(:, ch));
end
    
if length(ch) > 1
    if length(ch) == size(lfp.data, 2)
        % find lfp channel with highest SNR for events
        pow = fastrms(sig, 15);
        mEvent = mean(pow);
        meEvent = median(pow);
        mmEventRatio = mEvent ./ meEvent;
        mmEventRatio(mEvent < 1) = 0;
        mmEventRatio(meEvent < 1) = 0;
        [~, ch] = max(mmEventRatio);
    else
        % average signal across channels
        sig = mean(sig, 2);
    end
end

% rectify
if size(sig, 2) > 1
    sig = sig(:, ch) .^ 2;
else
    sig = sig .^ 2;
end
% moving average
sig = movmean(sig, winLength);
% standerdize
sig = (sig - mean(sig(inter))) / std(sig(inter));

% prepare noise
if ~isempty(noise)
    if ~isempty(passband)
        noise = filterLFP(double(lfp.data(:, noise)), 'type', type,...
            'passband', passband, 'order', order, 'graphics', false);
    else
        noise = double(lfp.data(:, noise));
    end
    noise = noise .^ 2;
    noise = movmean(noise, winLength);
    noise = (noise - mean(noise(inter))) / std(noise(inter));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUGGING; check signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% idx = round(60 * 60 * lfp.fs) : round(65 * 60 * lfp.fs);
% figure
% plot(lfp.data(idx))
% hold on
% plot(sig(idx));
% xlim([1 500000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detect events by thresholding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thresholded = sig > thr(1);
start = find(diff(thresholded) > 0);
stop = find(diff(thresholded) < 0);

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
nEvents = length(temp);

if isempty(temp)
	disp('Detection by thresholding failed');
	return
else
	disp(['After thresholding: ' num2str(nEvents) ' events.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge events if inter-event period is too short
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp2 = [];
event = temp(1, :);
for i = 2 : nEvents
	if temp(i, 1) - event(2) < interDur
		event = [event(1) temp(i, 2)];
	else
		temp2 = [temp2; event];
		event = temp(i, :);
	end
end
temp2 = [temp2; event];
nEvents = length(temp2);

if isempty(temp2)
	disp('Event merge failed');
	return
else
	disp(['After merging: ' num2str(nEvents) ' events.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discard events by duration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dur = temp2(:, 2) - temp2(:, 1);
idxMax = dur > maxDur * lfp.fs;
idxMin = dur < minDur * lfp.fs;
idx = idxMax | idxMin;
temp2(idx, :) = [];
nEvents = length(temp2);

if isempty(temp2)
	disp('duration test failed.');
	return
else
	disp(['After testing duration: ' num2str(nEvents) ' events.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exclude events that appear on noise channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bad = [];
if ~isempty(noise)           
	excluded = zeros(nEvents, 1);
	previous = 1;
	for i = 1 : nEvents
		if any(noise(temp2(i, 1) : temp2(i, 2)) > thr(2))
			excluded(i) = 1;
        end
	end
	bad = temp2(logical(excluded), 1);
	temp2 = temp2(~logical(excluded), :);
    nEvents = length(temp2);

    if isempty(temp2)
        disp('noisy channel test failed.');
        return
    else
        disp(['After noisy channel test: ' num2str(nEvents) ' events.']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exclude ripple-like events that are due to movement artifacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if emgThr   
    % load EMG file
    [~, basename, ~] = fileparts(basepath);
    if exist([basename '.emg.mat'])
        load([basename '.emg.mat']) 
    elseif exist([basename '.EMGfromLFP.mat'])
        load([basename '.EMGfromLFP.mat'])
    else
        warning('EMG file not found. Calculating EMG from LFP')
        emg = getEMGfromLFP(double(lfp.data(:, [1 4 7 10])), 'fs', lfp.fs, 'graphics', false, 'emgFs', 2);
    end
    
    % find artifacts
    excluded = zeros(nEvents, 1);
    for i = 1 : nEvents
        [~, idx] = min(abs(emg.timestamps - (temp2(i, 1) / lfp.fs)));
       if emg.data(idx) > emgThr
           excluded(i) = 1;           
       end
    end
    bad = sortrows([bad; temp2(logical(excluded), 1)]);
    temp2 = temp2(~logical(excluded), :);
    nEvents = length(temp2);
    
    if isempty(temp2)
        disp('EMG noise removal failed.');
        return
    else
        disp(['After EMG noise removal: ' num2str(nEvents) ' events.']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find peak power for each event
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

peakPower = zeros(nEvents, 1);
peakPos = zeros(nEvents, 1);
for i = 1 : nEvents
	[peakPower(i), peakPos(i)] = max(abs(sig([temp2(i, 1) : temp2(i, 2)])));
    peakPos(i) = peakPos(i) + temp2(i, 1) - 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange struct output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

events.preset = preset;
events.env = sig;
events.stamps = temp2;
events.peaks = peakPos;
events.power = peakPower;
events.noise = bad;
events.interval = interval;
events.ch = ch;
events.fs = lfp.fs;
% events.binary = zeros(size(lfp.timestamps));
% for i = 1 : nEvents
%     events.binary(temp2(i, 1) : temp2(i, 2)) = 1;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveVar   
    [~, filename] = fileparts(basepath);
    save([basepath, '\', filename, '.LFPevents.mat'], 'events')
end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    figure
    
    % standerd power and threshold
    subplot(2, 1, 1)
    plot(lfp.timestamps, sig);
    hold on
    axis tight
    Xlim = xlim;
    plot(Xlim, [thr(2) thr(2)])
    xlabel('Time [m]')
    ylabel('Standerdized power')
    box off
    set(gca, 'TickLength', [0 0])
    
    % example of detected ripple
    subplot(2, 1, 2)
    plot(lfp.timestamps, lfp.data(:, ch));
    Ylim = ylim;
    hold on
    for i = 1 : length(events.stamps)
        plot([events.stamps(i, 1) / lfp.fs, events.stamps(i, 1) / lfp.fs], Ylim, 'g')
        plot([events.peaks(i) / lfp.fs, events.peaks(i) / lfp.fs], Ylim, 'k')
        plot([events.stamps(i, 2) / lfp.fs, events.stamps(i, 2) / lfp.fs], Ylim, 'r')
    end
    [~, idx] = min(abs(events.power - mean(events.power)));
    idx = events.peaks(idx);
    axis tight
    xlim([idx / lfp.fs - 2, idx / lfp.fs + 2])
    xlabel('Time [s]')
    box off
    set(gca, 'TickLength', [0 0])
end

end

% EOF

