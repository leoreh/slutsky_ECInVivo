function [CaResults] = Ca_EventDetector(data,varargin)
% this funciton calc the dFF and finds events in the data. you can run it on any Ca2+ signal. it runs on each trace indevidually.  
% INPUT: structure with Ca2+ trace in a colm.  
% OUTPUT: 

pnames = {'params'};
dflts = {[]};
[params] = internal.stats.parseArgs(pnames,dflts,varargin{:});

%0) load params if not provided 
if isempty(params) 
%filtering params 
params.humWin = 10; 
params.smoothWin = 10; 
params.poly = 3; 
params.sgfwindow = 11; 
params.msadjustWINfactor = 0.1; % this is the relative size of the wind to the segmet size the user inputs later (0.1 = 10%)

%find peaks params 
params.method = 'halfprom';
params.minPeak = 0.12; % change this to play with events detection
params.minPeakDistance = 5; % change for mito 
params.minpeakProm = 0.02; 

else 
end

% 1) adjust the baseline for the whole trace  
timeVector = 1:1:length(data);
msadjustWIN =length(data)*params.msadjustWINfactor; 
close all 
dataAdjust = msbackadj(timeVector',data,'WindowSize',msadjustWIN,'ShowPlot',1); % use spline to adjust basline
answer = input('is the baseline ok? (type 1 to redo with diff win size, or enter to continue): ');
while answer == 1 
adjRatio = input('set new ratio (e.g 0.2): '); 
msadjustWIN =length(data)*adjRatio; 
close all
dataAdjust = msbackadj(timeVector',data,'WindowSize',msadjustWIN,'ShowPlot',1); % use spline to adjust basline
answer = input('is the baseline ok? (type 1 to redo with diff win size, or enter to continue): ');
end 
dataAdjust(dataAdjust<0) = 0; % turn all negative values to zero

% 2) ask user to choose frames to analyze
maxFrame = length(dataAdjust); 
close all 
plot(dataAdjust,'k'); set(gcf,'Color','White'); title('Choose Frames:')
frames = input(['look at the plot and choose frames to analyze from 1 to ', num2str(maxFrame), ' (e.g [1,120]): ']); 
dataSeg = dataAdjust(frames(1):frames(2)); 

% 3) Smooth the data in diff ways, calc. dF/F using percentile 
    %3.0) non filtered data is also a possiblity 
    FiltData(:,1) = dataSeg; 
    %3.1) Humming wnindow filter
    humFilt = hamming(params.humWin); % creates a vector for hamming window values. 
    normHumFilt= humFilt/sum(humFilt);	% Normalize window area under curve to 1
    FiltData(:,2) = filtfilt(normHumFilt,1,dataSeg); % filter cyto (filtfilt runs it forward and backwords so there is not phase shift)
    %3.2) smooth 
    FiltData(:,3) = smooth(dataSeg,params.smoothWin,'loess'); 
    %3.3) sg filter 
    FiltData(:,4) = sgolayfilt(dataSeg,params.poly,params.sgfwindow);
    %3.4) dFF
    params.prctile = input('set percentile for dF/F, recommended 20 for cyto, 5 for mito: ');     
    F = prctile(data(frames(1):frames(2)),params.prctile);
    dFF = FiltData/F; 

% 4) find events for each filter type and asks you to choose the best filter

close all;  
for filtType = 1:4
subplot(4,1,filtType)
time = 1:1:length(dFF(:,filtType));
findpeaks(dFF(:,filtType),time,'MinPeakHeight',params.minPeak,'MinPeakDistance',params.minPeakDistance,'MinPeakProminence',params.minpeakProm,'WidthReference',params.method,'Annotate','extents'); 
title(['filter number ',num2str(filtType)]); 
end

choosenFilt = input('the best is: 1.none, 2.humming, 3.smooth(loess), 4.sg ? (type number): ');

switch choosenFilt 
    case 1 
    [Events.pks, Events.locs, Events.width, Events.prom] = findpeaks(dFF(:,1),time,'MinPeakHeight',params.minPeak,'MinPeakDistance',params.minPeakDistance,'MinPeakProminence',params.minpeakProm,'Annotate','extents'); 
    CaResults.Filter = 'none';
    case 2 
    [Events.pks, Events.locs, Events.width, Events.prom] = findpeaks(dFF(:,2),time,'MinPeakHeight',params.minPeak,'MinPeakDistance',params.minPeakDistance,'MinPeakProminence',params.minpeakProm,'Annotate','extents');      
    CaResults.Filter = 'Humming';
    case 3 
    [Events.pks, Events.locs, Events.width, Events.prom] = findpeaks(dFF(:,3),time,'MinPeakHeight',params.minPeak,'MinPeakDistance',params.minPeakDistance,'MinPeakProminence',params.minpeakProm,'Annotate','extents');     
    CaResults.Filter = 'smooth_loess';
    case 4 
    [Events.pks, Events.locs, Events.width, Events.prom] = findpeaks(dFF(:,4),time,'MinPeakHeight',params.minPeak,'MinPeakDistance',params.minPeakDistance,'MinPeakProminence',params.minpeakProm,'Annotate','extents');       
    CaResults.Filter = 'SGF';
    otherwise 
        disp('ok try again with diff params for eventDetection!'); 
end 
        
        
% 5) now translate the events back to unfiltered data so the peaks are not cut off  
params.framesBack = 5; % here we choose how many frames back do i look for the peak. play with this for mito to work. 
unfiltered = dFF(:,1);


for eventNum = 1:size(Events.locs,2)   
backpk = unfiltered(Events.locs(eventNum)-params.framesBack:Events.locs(eventNum));   
locAdjust = (params.framesBack+1)-find(backpk == max(backpk)); % this calcualtes how many fatmes to go back to find the peak in the unfiltered data. 
CorrectEvents(eventNum,1) = Events.locs(eventNum)-locAdjust; 
CorrectEvents(eventNum,2) = unfiltered(CorrectEvents(eventNum,1));   
end 

%  a loop here to see if the postion correction is right. if not user he can cahnge params.framesback
close all; set(gcf,'Color','White')
plot(dFF(:,1),'k'); hold on; scatter(CorrectEvents(:,1),CorrectEvents(:,2),'r'); title('look if peaks are in the correct postion')
answer_framesback = input(['is the lag correction good? it now ',num2str(params.framesBack),' frames back. adjust it? (press y to change)'],'s');
while  answer_framesback == 'y'
     adjFramesBack = input('How many frames back to go? (defult was 5)');
     params.framesBack = adjFramesBack; 
for eventNum = 1:size(Events.locs,2)   
backpk = unfiltered(Events.locs(eventNum)-params.framesBack:Events.locs(eventNum));   
locAdjust = (params.framesBack+1)-find(backpk == max(backpk)); % this calcualtes how many fatmes to go back to find the peak in the unfiltered data. 
CorrectEvents(eventNum,1) = Events.locs(eventNum)-locAdjust; 
CorrectEvents(eventNum,2) = unfiltered(CorrectEvents(eventNum,1));   
end 
close all; set(gcf,'Color','White')
plot(dFF(:,1),'k'); hold on; scatter(CorrectEvents(:,1),CorrectEvents(:,2),'r'); title('look if peaks are in the correct postion')
answer_framesback = input(['is the lag correction good? it now ',num2str(params.framesBack),' frames back. adjust it? (y/n)'],'s');
end 


close all
%6)  go over each event and decide if you want keep it:
WIN = 15; % number of frame +/- to look at
i = 1;
selectedEvents = []; 

for eventNum = 1:size(CorrectEvents,1)
    
loc = CorrectEvents(eventNum,1); 
pk = CorrectEvents(eventNum,2);
STR = loc-WIN;
END = loc+WIN;

set(gcf,'Color','White')
subplot(2,1,1) % this just plots the whole trace
plot(dFF(:,1),'k'); hold on; scatter(loc,pk,'ro'); 
title('full trace'); 

% solving edge cases for locations 
if STR < 1 
 miniLoc = loc; 
    STR = 1;
else 
    miniLoc = WIN+1; 
end 

if END > length(dFF(:,1)) 
 END = length(dFF(:,1)); 
else 
end 
 
subplot(2,1,2) % this is the plot of the even to choose 
plot(dFF(STR:END,1),'k'); hold on; scatter(miniLoc,pk,'r*'); 
title(['Event number ',num2str(eventNum)]); 


eventStatus = input('save event? (1 = yes, 2 = no): '); 
 if eventStatus == 1  
  selectedEvents(i,1) = loc;
  selectedEvents(i,2) = pk;
  i = i+1; 
 else 
 end
 close all 
end    

% now give a report back 
disp(['analysis finshed (see final plot).', num2str(length(selectedEvents)),' events approved. CaResults holds the data and params']); 
set(gcf,'Color','White')
plot(dFF(:,1),'k'); hold on; scatter(selectedEvents(:,1),selectedEvents(:,2),'r'); 
title('Final Result. All the selected events (in red)')

%7) now arrange the data in one structure 
CaResults.dFF = dFF; % the dFF traces we calcaulted together
CaResults.params = params; % saves all the params you used for later ref
CaResults.Events = CorrectEvents; % all events projected on unfiltered data 
CaResults.Selected_Events = selectedEvents; % only the selected events 
CaResults.frames = frames; % the frames we analzyed togehter 

    
end % final end 
    

