function bs = findBS(varargin)

% detects ripples from LFP signal. based on bz_FindRipples. main
% differences are (1) moving average implemented via movmean, (2) detects
% negative peak after elimination, (3) selects best channel automatically.
% more information is availble in PreProcessing.docx
% 
% INPUT
%   lfp         struct (see getLFP.m)
%   ch          channel to detect ripples. if empty will select by RMS
%   noise       channel to remove artifacts {[]}
%   emgThr      threshold for EMG exclusion {0.9}
%   interval    interval used for standadization {[0 Inf]}
%   basepath    recording session path {pwd}
%   graphics    plot figure {1}.
%   saveVar     save variable {1}.
% 
% OUTPUT
%   ripples         struct with fields:
% 
% CALLS
%   fastrms
% 
% TO DO LIST:
%   compare with original script on various data
%   track and save exclusions
%   add raw lfp to graphics
% 
% 01 may 19 LH. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'lfp', @isstruct)
addParameter(p, 'ch', [], @isnumeric)
addParameter(p, 'noise', [], @isnumeric)
addParameter(p, 'emgThr', 0.9, @isnumeric);
addParameter(p, 'interval', [0 Inf], @isnumeric)
addParameter(p, 'basepath', pwd, @isstr);
addParameter(p, 'graphics', true, @islogical)
addParameter(p, 'saveVar', false, @islogical);

parse(p, varargin{:})
lfp = p.Results.lfp;
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

% prepare output
bs = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
passband = [2 100];           % [Hz]
order = 3; 
type = 'butter';
winLength = 11;                 % [samples]    

interRipple = 0.2 * lfp.fs;    % minimum between ripple events [s]
thr = [0.5 1];                    % low and high thresholding value [z-score]
maxDur = 4000;                   % maximum ripple duration [ms]
minDur = 0.0001;                  % minimum ripple duration [ms]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filter    
sigfilt = filterLFP(double(lfp.data), 'type', type,...
    'passband', passband, 'order', order, 'graphics', false);
    
% find lfp channel with highest SNR for ripples
if isempty(ch)    
    pow = fastrms(sigfilt, 15);
    mRipple = mean(pow);
    meRipple = median(pow);
    mmRippleRatio = mRipple ./ meRipple;   
    mmRippleRatio(mRipple < 1) = 0;
    mmRippleRatio(meRipple < 1) = 0;   
    [~, ch] = max(mmRippleRatio);
end

% rectify
signorm = sigfilt(:, ch) .^ 2;
% moving average
signorm = movmean(signorm, winLength);
% standerdize
signorm = (signorm - mean(signorm(inter))) / std(signorm(inter));

% prepare noise
if ~isempty(noise)
    noise = sigfilt(:, noise) .^ 2;
    noise = movmean(noise, winLength);
    noise = (noise - mean(noise(inter))) / std(noise(inter));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detect ripples by thresholding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thresholded = signorm > thr(1);
start = find(diff(thresholded) > 0);
stop = find(diff(thresholded) < 0);

% exclude last ripple if it is incomplete
if length(stop) == length(start) - 1
	start = start(1 : end - 1);
end
% exclude first ripple if it is incomplete
if length(stop) - 1 == length(start)
    stop = stop(2 : end);
end
% correct special case when both first and last ripples are incomplete
if start(1) > stop(1)
	stop(1) = [];
	start(end) = [];
end

firstPass = [start, stop];
nRippels = length(firstPass);

if isempty(firstPass)
	disp('Detection by thresholding failed');
	return
else
	disp(['After detection by thresholding: ' num2str(nRippels) ' events.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge ripples if inter-ripple period is too short
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% option 1
% idx = find(firstPass(2 : end, 1) - firstPass(1 : end - 1, 2) >= interRipple);
% secondPass = [firstPass(1, :); firstPass(idx + 1, 1), firstPass(idx + 1, 2)];

%%% option 2 (original)
secondPass = [];
ripple = firstPass(1, :);
for i = 2 : nRippels
	if firstPass(i,1) - ripple(2) < interRipple
		ripple = [ripple(1) firstPass(i, 2)];
	else
		secondPass = [secondPass; ripple];
		ripple = firstPass(i, :);
	end
end
secondPass = [secondPass; ripple];
nRippels = length(secondPass);

if isempty(secondPass)
	disp('Ripple merge failed');
	return
else
	disp(['After ripple merge: ' num2str(nRippels) ' events.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discard ripples by duration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thirdPass = secondPass;
dur = thirdPass(:, 2) - thirdPass(:, 1);
idxMax = dur > maxDur * lfp.fs;
idxMin = dur < minDur * lfp.fs;
idx = idxMax | idxMin;
thirdPass(idx, :) = [];
nRippels = length(thirdPass);

if isempty(thirdPass)
	disp('duration test failed.');
	return
else
	disp(['After duration test: ' num2str(nRippels) ' events.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exclude ripple-like events that appear on noise channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bad = [];
if ~isempty(noise)           
	excluded = zeros(nRippels, 1);
	previous = 1;
	for i = 1 : nRippels
		if any(noise(thirdPass(i, 1) : thirdPass(i, 2)) > thr(2))
			excluded(i) = 1;
        end
	end
	bad = thirdPass(logical(excluded), 1);
	thirdPass = thirdPass(~logical(excluded), :);
    nRippels = length(thirdPass);

    if isempty(thirdPass)
        disp('noisy channel test failed.');
        return
    else
        disp(['After noisy channel test: ' num2str(nRippels) ' events.']);
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
    excluded = zeros(nRippels, 1);
    for i = 1 : nRippels
        [~, idx] = min(abs(emg.timestamps - (thirdPass(i, 1) / lfp.fs)));
       if emg.data(idx) > emgThr
           excluded(i) = 1;           
       end
    end
    bad = sortrows([bad; thirdPass(logical(excluded), 1)]);
    thirdPass = thirdPass(~logical(excluded), :);
    nRippels = length(thirdPass);
    
    if isempty(thirdPass)
        disp('EMG noise removal failed.');
        return
    else
        disp(['After EMG noise removal: ' num2str(nRippels) ' events.']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discard ripples with low peak power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thirdPass = [];
peakPower = [];
for i = 1 : nRippels
	[val] = max(signorm([secondPass(i, 1) : secondPass(i, 2)]));
	if val > thr(2)
		thirdPass = [thirdPass; secondPass(i, :)];
		peakPower = [peakPower; val];
    end
end
nRippels = length(thirdPass);

if isempty(thirdPass)
	disp('Peak thresholding failed.');
	return
else
	disp(['After peak thresholding: ' num2str(nRippels) ' events.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detect negative peak position for each ripple
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

peakPos = zeros(nRippels, 1);
for i = 1 : nRippels
	[~, idx] = min(sigfilt(thirdPass(i, 1) : thirdPass(i, 2)));
	peakPos(i) = idx + thirdPass(i, 1) - 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange struct output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bs.timestamps = [lfp.timestamps(thirdPass(:, 1)), lfp.timestamps(thirdPass(:, 2))];
bs.peaks = lfp.timestamps(peakPos);
bs.power = peakPower;
bs.noise = bad;
bs.interval = interval;
bs.ch = ch;
bs.fs = lfp.fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveVar   
    [~, filename] = fileparts(basepath);
    save([basepath, '\', filename, '.ripples.mat'], 'ripples')
end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    figure
    
    % standerd power and threshold
    subplot(2, 1, 1)
    plot(lfp.timestamps / 60, signorm);
    hold on
    axis tight
%     ylim([0 thr(2) * 3]);
    Xlim = xlim;
    plot(Xlim, [thr(2) thr(2)])
    xlabel('Time [m]')
    ylabel('Standerdized ripple power')
    box off
    set(gca, 'TickLength', [0 0])
    
    % example of detected ripple
    subplot(2, 1, 2)
    plot(lfp.timestamps, lfp.data(:, ch));
%     ylim([-200 200])
    Ylim = ylim;
    hold on
    for i = 1 : length(bs.timestamps)
        plot([bs.timestamps(i, 1) bs.timestamps(i, 1)], Ylim, 'g')
        plot([bs.peaks(i) bs.peaks(i)], Ylim, 'k')
        plot([bs.timestamps(i, 2) bs.timestamps(i, 2)], Ylim, 'r')
    end
    [~, x] = max(bs.power);
    x = bs.timestamps(x, 2);
%     xlim([x - 0.25, x + 0.25])
    xlabel('Time [s]')
    box off
    set(gca, 'TickLength', [0 0])   
end

end

% EOF

