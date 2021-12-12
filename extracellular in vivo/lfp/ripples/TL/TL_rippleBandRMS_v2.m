function [beztRipChans] = TL_rippleBandRMS_v2(pathname , saveRMS)

% v2: loads expInfo instead of using excel spreadsheet. DOES NOT WORK FOR
% MULTIPLE SETS YET

%% Find which recordings are baseline

if ~strcmp(pathname(end) , filesep)
    pathname = [pathname , filesep];
end

% Find recording folders;
% [recFols] : cell array of recording folders
d = dir(pathname);
d = {d([d.isdir]).name};
i = cell2mat(cellfun(@(z) ~strcmp(z(1) , '.') , d , 'uniformoutput' , false));
recFols = [d(i)];
clear d i;

bL = nan(length(recFols) , 1);
for r = 1 : length(recFols)
    
    recPath = [pathname , recFols{r} , filesep];
    
    % Load excel file with experiment info
    d = dir([recPath , '*expInfo.mat']);
    try
    fname = d.name;
    load([recPath , fname]);
    catch
        bL(r) = NaN;
    end
    bL(r) = ~isempty(strfind(expInfo.notes{1} , 'aseline'));    
end

baselineFols = [recFols(logical(bL))]; %%%%%%%%%%%%%%%%%%%%%%
% exc = [exc(logical(bL))];
clear bL recFols;

%% Calculate rms from each channel across 6 hour epochs during baseline, and find the best

RMS = [];

for b = 1 : length(baselineFols)
    currentPath = [pathname , baselineFols{b} , filesep];
    
    % Check directory in pathname ---------------------------------------------
    % [fname] : names of files in pathname, including their suffix (eg .mat ,.xls)
    temp = dir(currentPath);
    filenames = {temp(~[temp.isdir]).name}';
    clear temp;
    
    % Load sleep states and convert to conformity
    sleep = TL_extractSleep(currentPath)
    
    % Calculate ranges for each chunk of lfp extraction
    ranges = 0 : 60*60 : length(sleep);
    
    ch = expInfo.channels{1};
    
    for r = 1 : length(ranges) - 1
        currRow = size(RMS , 1) + 1;
        for c = 1 : length(ch)
            lfp = TL_getLFP('basepath' , currentPath , 'ch' , ch(c) , 'fs' , 1250 , ...
                'extension' , 'lfp' , 'savevar' , true , 'forceL' , true , 'interval' , [ranges(r) , ranges(r+1)]);
            
            % Run butterworth Filter
            [c1 c2] = butter(5,[100 250]/(1250/2),'bandpass');
            filtD = filtfilt(c1,c2,double(lfp.data));
            
            % Extract Only NREM states from the filtered Data
            % this is not exactly current because of the concatenations
            % between different NREM states but all channels should have
            % the same number of concatenations at least
            rEnd = ranges(r+1);
            if r == length(ranges) - 1 & rEnd > length(sleep)
                rEnd = length(sleep) - 1;
            end
            currSleep = sleep([ranges(r) :  rEnd] + 1);
            currSleep = currSleep == 4;
            currSleep = repelem(currSleep , 1250);
            if length(currSleep) > length(lfp.data)
            currSleep = currSleep(1:length(lfp.data));
            end
            d = find(diff(currSleep) == 1) + 1;
            if size(d,1) == 1,
                d = d';
            end
            if currSleep(1) == 1
                d = [1 ; d];
            end
            dd = find(diff(currSleep) == - 1) - 1;
            if size(dd,1) == 1,
                dd = dd';
            end
            if currSleep(end) == 1
                dd = [dd ; length(currSleep)];
            end
            epochs = [d , dd];
            clear d dd;
            tempRMS = [];
            for e = 1 : size(epochs , 1)
                tempRMS = [tempRMS , repelem(rms(filtD(epochs(e,1):epochs(e,2))) , diff(epochs(e,:)))];
            end
            RMS(currRow,c) = mean(tempRMS); %rms(filtD(logical(currSleep)));
            clear tempRMS;
        end
        %%
    end
end

% Now find the best channel
[~ , bezt] = sort(mean(RMS,1),'descend');
beztRipChans.RMS = RMS;
beztRipChans.rank = bezt;

% need to find which tetrode this best channel is on
% ch = expInfo.channels{1};
% rm = expInfo.channels{1};
% for r = 1 : length(rm)
%     i = ch>=rm(r)
%     ch(i) = ch(i) + 1;
% end
% ch = sort([ch,rm],'ascend');
b = bezt + sum(bezt >= expInfo.channels{1});
tet = floor(b/4);
beztRipChans.tet = tet;

if saveRMS
    savePath = save([pathname , 'beztRipChans'],'beztRipChans');
end