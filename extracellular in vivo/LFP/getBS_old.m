function bs = getBS(varargin)

% detects burst events from LFP using the std of the signal.
% 
% INPUT
%   sig         signal for detection
%   fs          sampling frequency
%   thr         thr for std separation. if empty then determined via
%               bimodal distribution of data
%   basepath    recording session path {pwd}
%   graphics    logical. plot figure {1}.
%   saveVar     logical. save variable {1}.
% 
% OUTPUT
%   events         struct with fields:
%       stamps      n x 2 mat where n is the number of events and the
%                   columns represent the start and end of each event [samps]
%       peaks       the timestamp for each peak [samps] 
%       peakPower   the voltage at each peak [uV]
%       <params>    as in input
% 
% CALLS
%   remNear
% 
% 22 nov 19 LH. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'sig', [], @isnumeric)
addParameter(p, 'fs', 1250, @isnumeric)
addParameter(p, 'thr', [], @isnumeric)
addParameter(p, 'basepath', pwd, @isstr);
addParameter(p, 'graphics', true, @islogical)
addParameter(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
sig = p.Results.sig;
fs = p.Results.fs;
thr = p.Results.thr;
basepath = p.Results.basepath;
graphics = p.Results.graphics;
saveVar = p.Results.saveVar;

% params
winmax = 0.8 * fs;                  % window for moving max        
minDur = 1;                         % minimum event duration [s]
maxDur = 300;                        % maximum event duration [s]
stdwin = round(minDur * 2 * fs);    % bin width for std calculation
dist = 0.2 * fs;                    % distance for removing adjacent peaks
interDur = 2 * fs;                  %  minimum time between events [samples]
        
% prepare output
bs = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep signal and detect events
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filter
% x = filterLFP(sig, 'fs', fs, 'passband', [0 20], 'order', 6,...
%     'type', 'butter', 'dataOnly', true, 'graphics', false,...
%     'saveVar', false);

% threshold for bs separation is determined according to the bimodel
% distribution of std values
x = movstd(double(sig), stdwin);
if isempty(thr)
    thr = sepBimodel('x', x, 'lognorm', true, 'graphics', true);
end

% detection
b = x > thr;

% find start
start = find(diff(b) > 0);
% start = movmax(x, [winmax, 0]);
% start = find(diff(start) > thr)';
% start = remNear('x', start, 'dist', dist, 'flip', false);

% find end
stop = find(diff(b) < 0);
% stop = movmax(x, [0, winmax]);
% stop = find(diff(stop) < thr)';
% stop = remNear('x', stop, 'dist', dist, 'flip', true);

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

iei = temp(2 : end, 1) - temp(1 : end - 1, 2);      % inter-event interval
while min(iei) < interDur
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
    temp = temp2;
    iei = temp(2 : end, 1) - temp(1 : end - 1, 2);
end

nEvents = length(temp);

if isempty(temp)
	disp('Event merge failed');
	return
else
	disp(['After merging: ' num2str(nEvents) ' events.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discard events by duration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dur = temp(:, 2) - temp(:, 1);
idxMax = dur > maxDur * fs;
idxMin = dur < minDur * fs;
idx = idxMax | idxMin;
temp(idx, :) = [];
nEvents = length(temp);

if isempty(temp)
	disp('duration test failed.');
	return
else
	disp(['After testing duration: ' num2str(nEvents) ' events.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find peak power for each event
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

peakPower = zeros(nEvents, 1);
peakPos = zeros(nEvents, 1);
for i = 1 : nEvents
	[peakPower(i), peakPos(i)] = max(abs(sig([temp(i, 1) : temp(i, 2)])));
    peakPos(i) = peakPos(i) + temp(i, 1) - 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange struct output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bs.stamps = temp;
bs.peaks = peakPos;
bs.power = peakPower;
bs.fs = fs;
bs.thr = thr;
% bs.binary = zeros(length(b), 1);
% for i = 1 : size(bs.stamps, 1)
%     bs.binary(bs.stamps(i, 1) : bs.stamps(i, 2)) = 1;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveVar   
    [~, filename] = fileparts(basepath);
    save([basepath, '\', filename, '.bs.mat'], 'bs')
end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    figure
    
    % standerd power and threshold
    subplot(2, 1, 1)
    plot(lfp.timestamps, x);
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

