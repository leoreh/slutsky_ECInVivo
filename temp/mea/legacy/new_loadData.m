%export .plx file with: channel, unit, timestamp, PC1, PC2, peak, valley, raw waveform data in microvolts

 

function data = new_loadData(sessionKey, bSilent)
  
    persistent local_data
    persistent local_session_key    
   
   
    
    if ~isempty(local_data) && strcmp(sessionKey, ...
				      local_session_key)==1
      data = local_data;
      return
    end
      
    if ~exist('bSilent', 'var')
        bSilent = false;
    end    
    sessionConfig = getSessionConfig(sessionKey);
    [dataFile,dataDir] = getSessionDataPath(sessionKey);
    dataFile=dataFile.filePath;
    %cd(dataDir); 
    if ~bSilent
        fprintf('Loading %s (from %s)\n', sessionKey, dataFile)
    end
    data = struct;
    data.sessionKey = sessionKey;
    data.dataFile = dataFile;
    data.nBaselineHours = sessionConfig.nBaselineHours;
    fileData = load(dataFile);
    data = addFromFile(data,fileData);
    
    
    local_data = data; 
    local_session_key = sessionKey;
    
    
end



% ==========================================
function data = addFromFile(data, fileData)% fileData is the m.file created form plexon
    shortUnitThreshold = 0.09; % ignore units that don't have spikes for at least this number of hours
    
    data.channelNumbers = [];
    data.channels = cell(1,0);
    data.unitIds = {};
    data.unitSpikeTimes = cell(1,0);
    data.unitFullHours = [];
    data.unitPC = [];
    data.unitFreq = zeros(120); %%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<changed from 120
    %data.unitFreq = data.unitFreq+99;
    data.normUnitFreq = zeros(120);%%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<changed from 120
    %data.normUnitFreq = data.normUnitFreq +99;
    data.allUnitPC = [];
    data.unitPeakValley = [];
    data.allunitPeakValley =[];
    nChannels = 0;
    nUnits = 0;
    nShortUnits = 0;
    nUnstableUnits = 0;
    UnstableUnits_unit_num = [];
    UnstableUnits_channel = []; 
    nSlowUnits = 0;
    SlowUnits_unit_num = [];
    SlowUnits_channel = [];
    allFields = fieldnames(fileData);
    SegTime = 1200;%%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<changed from 1200
  
    for iField = 1:length(allFields)
        f = allFields{iField};
        channel_data = fileData.(f);
        if isempty(channel_data) % if channel is empty dont go further, go back to for and next channel
            continue
        end
        nChannels = nChannels + 1;
        chNum = str2double(f(4:end));        
        data.channelNumbers(nChannels) = chNum;
        data.channels{nChannels} = channel_data;
        units = unique(channel_data(:,2));
        for iUnit = 1:length(units)
            unitNum = units(iUnit);
            unitRows = channel_data(:,2) == unitNum;
            unitTimes = channel_data(unitRows,3);
            nFullHours = round(max(unitTimes)/SegTime); % "hour" is only 20 minutes. Count the hour if it had a spike its 2nd half
            bins=(0:SegTime:nFullHours*SegTime);
            unitFreq=histc(unitTimes,bins);
            unitFreq=reshape(unitFreq,[],1);
            unitFreq=unitFreq(1:nFullHours,:)/SegTime;
%             normUnitFreq=100*unitFreq/mean(unitFreq(1:data.nBaselineHours));
           % distFromZero= sqrt((mean(channel_data(unitRows,4)))^2 + mean(channel_data(unitRows,5))^2);
           % unitPC= [mean(channel_data(unitRows,4)), mean(channel_data(unitRows,5)), distFromZero];
           % unitPeakValleyDiff= max(mean(channel_data(unitRows,8:37))) - min(mean(channel_data(unitRows,8:37)));
           % unitPeakValley = [ max(mean(channel_data(unitRows,8:37))),  min(mean(channel_data(unitRows,8:37))), unitPeakValleyDiff];
%             if  nFullHours < shortUnitThreshold || numel(find(unitFreq(12:numel(unitFreq),:) == 0))>6 ...
%                     || numel(find(unitFreq(5:numel(unitFreq),:)))<10;
% 
%                 nShortUnits = nShortUnits + 1;
%                  continue
%             end
%             if  isStable(unitTimes, data.nBaselineHours,0.1)% only units that changed a 1000 from 1 to 3 hour are unstable...?!?
%                 nUnstableUnits = nUnstableUnits + 1;
%                 UnstableUnits_unit_num(end+1) = unitNum;
%                 UnstableUnits_channel(end+1) = chNum;
%                 continue
%             end
          if  isSlow(unitTimes, data.nBaselineHours, 0.1) % note that i changed it from 0.1 to 0.01!!!!
                 nSlowUnits = nSlowUnits + 1;
                 SlowUnits_unit_num(end+1) = unitNum;
                 SlowUnits_channel(end+1) = chNum;
                 continue
            end
            nUnits = nUnits + 1;
            data.unitIds(nUnits,1:2) = {f unitNum};
            data.unitSpikeTimes{nUnits} = unitTimes;
%             data.unitPC{nUnits} = unitPC;
%             data.allUnitPC(nUnits,1:3) = unitPC;
%             data.unitPeakValley{nUnits} = unitPeakValley;
%             data.allunitPeakValley(nUnits,1:3) = unitPeakValley;
            data.unitFreq(1:nFullHours,nUnits) = unitFreq;
            %data.normUnitFreq(1:nFullHours,nUnits) = normUnitFreq;
           % data.unitWaveform{nUnits} = mean(channel_data(unitRows,8:37));

            data.unitFullHours(nUnits) = nFullHours;
        end        
    end
    data.networkFullHours = min(data.unitFullHours);
    data.maxUnitFullHours = max(data.unitFullHours);
    data.nChannels = nChannels;
    data.nUnits = nUnits;
    data.nDiscardedShortUnits = nShortUnits;
    data.nDiscardedUnstableUnits = nUnstableUnits;
    data.UnstableUnits_unit_num = UnstableUnits_unit_num;    
    data.UnstableUnits_channel = UnstableUnits_channel;    
    data.nSlowUnits = nSlowUnits;
    data.SlowUnits_unit_num = SlowUnits_unit_num;
    data.SlowUnits_channel = SlowUnits_channel;
    
end

% =========================================================
function bStable = isStable(spikeTimes, nBaselineHours, threshold)
     
ISIs = calcISIs(spikeTimes);
    fRate = @(x) 1/mean(x);
    rateStart = fRate(ISIs{1});
    rateEnd = fRate(ISIs{nBaselineHours});
    deltaRate = rateEnd-rateStart;
    change = deltaRate/rateStart;
    bStable = abs(change) < threshold;
end


% =========================================================
function bSlow = isSlow(spikeTimes, nBaselineHours, threshold)
    ISIs = calcISIs(spikeTimes);
    fRate = @(x) 1/mean(x);
    rateStart = fRate(ISIs{1});
    bSlow = rateStart < threshold;
end



