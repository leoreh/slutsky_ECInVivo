function CytoMitoEvents = Ca_EventCompare(dFF,CaEvents)
% this function  assumes that cyto preceds mito events. so looks for time locked events 
% we go mito event by mito event, and seek the correct cyto event. 
%INPUT1: dFF has colm 1 that is dFF of Cyto, and colm 2 that is dFF of mito (take them from CaResutls.dFF(:,1))
%INPUT2: CaEvents is a strucutre that has Events.cyto and Events.mito. each one has a list of all locs and pks in colms 1 and 2. take this from the CaResults.selectedEvents

winStart = 80; % play with this window sizes untill happy  
winEnd = 200;
i = 1; 
CytoMitoEvents = []; 

% perpare
Events.Mlocs = CaEvents.mito(:,1);
Events.Mpks =  CaEvents.mito(:,2);
Events.Clocs = CaEvents.cyto(:,1);
Events.Cpks =  CaEvents.cyto(:,2);
CaData.mito = dFF(:,2);
CaData.cyto = dFF(:,1);

for eventNum = 1:length(Events.Mlocs)
    Mloc = Events.Mlocs(eventNum);
    Cloc = Events.Clocs(find(Events.Clocs<Mloc,1,'last')); % find the closest preceding cyto event to the mito 
    CeventNum = find(Events.Clocs<Mloc,1,'last');
   
   figure('Position', [50 50 1000 550]) % this may need to chage according to screen
   subplot(2,1,1) 
   plot(CaData.mito)
   hold on; scatter(Mloc,Events.Mpks(eventNum),'b*')
   hold on; plot(CaData.cyto)
   hold on;  scatter(Cloc,Events.Cpks(CeventNum),'r*')
   
   subplot(2,1,2)
   plot(Mloc-winStart:Mloc+winEnd,CaData.mito(Mloc-winStart:Mloc+winEnd))
   hold on; scatter(Mloc,Events.Mpks(eventNum),'r*')
   hold on; plot(Mloc-winStart:Mloc+winEnd,CaData.cyto(Mloc-winStart:Mloc+winEnd))
   hold on;  scatter(Cloc,Events.Cpks(CeventNum),'b*')
   legend('mito','mito peak','cyto','cyto peak')
   
MaximAnswer = input('Save event? (y/n/e)  ','s'); 
if MaximAnswer == 'y'    
   CytoMitoEvents(i,1) = CeventNum; % save index for Cyto event 
   CytoMitoEvents(i,2) = eventNum; % save index for Mito event that corresponds 
   i = i+1; 
   close all  
elseif MaximAnswer == 'e' 
    disp('goodbye!')
    close all
    break 
else 
    disp('ok, on to the next'); 
    close all 
end 
clf 
end 
end 
    

