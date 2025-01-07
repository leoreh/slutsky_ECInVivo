function [cat] = TL_catMouse_v2(pathname , logIndx)

% INPUTS:
% [pathname] : string to subject folder (contains all date folders)
% [logIndex] : logical index of folders to use (in numerical order)


% Concatenates data from each day of recording for a given mouse, places in
% folder named "cat" . Must be run after TL_spikeWrapper & TL_rippleWrapper

% adds Delta/Theta ratio

% currently not appreciate multiple blocks in the same set (ignores gap)

temp = dir(pathname);
folnames = {temp([temp.isdir]).name};
i = cell2mat(cellfun(@(z) ~strcmp(z(1),'.') , folnames , 'uniformoutput' , false));
folnames = [folnames(i)];

if ~isempty(logIndx)
    folnames = [folnames(logical(logIndx))];
end

% TL not working for multiple sets of ripple files in the same folder

cat.MuaZbyDay = [];
cat.MuaHz = [];
cat.start = [];
cat.stop = [];
cat.ratio = [];
cat.spt = [];
cat.dayT = [];
cat.sleep = [];
cat.ripHz = [];
cat.MuaNonRipHz = [];
cat.ripSpkDayT = [];
cat.theta = [];
cat.delta = [];
cat.boutDayT = [];
cat.boutLabels = [];
for f = 1 : length(folnames)
    
    temp = dir([pathname , filesep , folnames{f} , filesep , 'ripples' , '*.mat']);
    load([pathname , filesep , folnames{f} , filesep , temp.name]);
    clear temp;
    
    load([pathname , filesep , folnames{f} , filesep , 'sleep.mat']);
    
    temp = dir([pathname , filesep , folnames{f} , filesep , '*spktimes.mat']);
    load([pathname , filesep , folnames{f} , filesep , temp.name]);
    clear temp;
    
    temp = dir([pathname , filesep , folnames{f} , filesep , '*Raw' , '*expInfo.mat']);
    load([pathname , filesep , folnames{f} , filesep , temp.name]);
    clear temp;
    
    
    % to find channels, pick the tetrode with the mose spikes, then pick
    % the channel with the most spikes on that tetrode
    [~ , tet] = max(cell2mat(cellfun(@(z) length(z) , spktimes , 'uniformoutput' , false)));
    ch = mode(spkch{tet});
    
    if strcmp(folnames{f} , '210723')
        [lfp , ~] = TL_getLFPtest('basepath' , [pathname , filesep , folnames{f} , filesep] , ...
            'type' , 'filt' , 'ch' , ch , 'fs' , 1250 , 'nchans' , length(expInfo.channels{1}) , ...
            'extension' , 'lfp' , 'savevar' , false , 'forceL' , true , 'interval' , [0 , 22*60*60]);
    else
        [lfp , ~] = TL_getLFPtest('basepath' , [pathname , filesep , folnames{f} , filesep] , ...
            'type' , 'filt' , 'ch' , ch , 'fs' , 1250 , 'nchans' , length(expInfo.channels{1}) , ...
            'extension' , 'lfp' , 'savevar' , false , 'forceL' , true);
    end
    
    [nremI , wakeI , theta , delta , ratio] = TL_deltaTheta(lfp , 1250 , [] , []);
    
    % only zscore with non"Bin" bouts
    labels = labels(1:length(ratio));
    
    setnum = find(cell2mat(cellfun(@(z) ~isempty(z) , expInfo.startTime , 'uniformoutput' , false)));
    setI = find(cell2mat(cellfun(@(z) z == setnum , expInfo.set , 'uniformoutput' , false)));
    
    
    % Extract spike times and convert to z score and raw time relative to start of exp
    m = cell2mat(cellfun(@(z) length(z) , spktimes , 'uniformoutput' , false));
    [~ , tet] = max(m);
    clear m;
    sptDays = (spktimes{tet}/24414.06)/(24*60*60) + expInfo.startTime{setI};
    
    % convert spt to zscore mua per second
    e = ([0 : 1 : expInfo.nsec{setI}])/(24*60*60) + expInfo.startTime{setI};
    
       muaBinnedHz = histcounts(sptDays , e);
    smMuaBinnedHz = smoothdata(muaBinnedHz , 'gaussian' , 5);
    
    if length(ratio) == length(smMuaBinnedHz)
        m = mean(smMuaBinnedHz);
        sd = std(smMuaBinnedHz);
        MuaZbyDay = (smMuaBinnedHz - m)/sd;
        DayT = e(1:end-1);
        start = expInfo.startTime{setI};
        stop = expInfo.startTime{setI} + expInfo.nsec{setI}/(24*60*60);
    else
        m = mean(smMuaBinnedHz(1:length(ratio)));
        sd = std(smMuaBinnedHz(1:length(ratio)));
        MuaZbyDay = (smMuaBinnedHz(1:length(ratio)) - m)/sd;
        muaBinnedHz = muaBinnedHz(1:length(ratio));
        DayT = e(1:length(ratio));
        start = expInfo.startTime{setI};
        stop = expInfo.startTime{setI} + length(MuaZbyDay)/(24*60*60);
    end

    sptDays = sptDays(sptDays < (spktimes{tet}/24414.06 + length(MuaZbyDay))/(24*60*60) + expInfo.startTime{setI});
    ripples.detect.lfpIndx = (ripples.detect.lfpIndx/1250)/(24*60*60) + start; % convert to true time
    keepRipI = ripples.detect.lfpIndx(:,2) < DayT(end);
          
    ripples.detect.DayT = ripples.detect.lfpIndx;
    ripples.detect = rmfield(ripples.detect,'lfpIndx');
    
                % Convert all ripple spiketimes to DayT
    ripSpkDayT = [];
    ripSpkCnt = [];
    for r = 1 : size(ripples.detect.DayT , 1)
        spkDayT = ripples.detect.spkT{r}/(24*60*60);
        if ~isempty(spkDayT)
            durDayT = ripples.detect.DayT(r,2) - ripples.detect.DayT(r,1);
            i = spkDayT > 0 & spkDayT < durDayT;
            ripSpkCnt(r) = sum(i);
            keepSpkT = spkDayT(i);
            ripSpkDayT = [ripSpkDayT; keepSpkT + ripples.detect.DayT(r,1)];
        end
    end
    
    ripBinnedHz = histcounts(ripSpkDayT , e);
    ripBinnedHz = ripBinnedHz(1:length(ratio));
    
    % Extract # non ripple spikes per seconnd
    MuaNonRipHz = muaBinnedHz - ripBinnedHz;
    
    fn = fieldnames(ripples.detect);
    
    d = round(24*60*60*(DayT - DayT(1)),0);
    r = round(floor(24*60*60*(ripSpkDayT - DayT(1))),0);
    for rr = 1 : length(r)
        fd = d == r(rr);
        ripSleep(rr) = labels(fd);
    end
    
    if f == 1
        catRip = ripples.detect;
        catRip.spkCnt = ripSpkCnt(keepRipI)';
    else
        catRip.spkCnt = [catRip.spkCnt; ripSpkCnt(keepRipI)'];
        for ff = 1 : length(fn)
            catRip.(fn{ff}) = [catRip.(fn{ff}) ; ripples.detect.(fn{ff})(keepRipI,:)];
            
        end
    end
    % Save bouts
    % TL this section also not appreciate timing for multiple blocks in the same set
    boutn = [expInfo.boutNames(setI)];
    boutt = [expInfo.EpochStartSec(setI)];
    boutlabels = {};
    boutDayT = [];
    for e = 1 : length(boutn)
        if iscell(boutn{e})
            for ee = 1 : length(boutn{e})
                boutlabels{end+1} = boutn{e}{ee};
                boutDayT(end+1) = start + boutt{e}(ee)/(24*60*60);
            end
        else
            boutlabels{end+1} = boutn{e};
            boutDayT(end+1) = start + boutt{e}/(24*60*60);
        end
        
    end
    
    cat.dayT = [cat.dayT, DayT];
    cat.MuaZbyDay = [cat.MuaZbyDay, MuaZbyDay;];
    cat.start = [cat.start, start];
    cat.stop = [cat.stop, stop];
    cat.theta = [cat.theta , theta];
    cat.delta = [cat.delta , delta];
    cat.spt = [cat.spt , sptDays'];
    cat.sleep = [cat.sleep , labels(1:length(ratio))'];
    cat.MuaHz = [cat.MuaHz , muaBinnedHz];
    cat.MuaNonRipHz = [cat.MuaNonRipHz , MuaNonRipHz];
    cat.ripHz = [cat.ripHz , ripBinnedHz];
    cat.ripSpkDayT = [cat.ripSpkDayT ; ripSpkDayT];
    cat.boutDayT = [cat.boutDayT , boutDayT];
cat.boutLabels = [cat.boutLabels , boutlabels];
    clear ripHz MuaNonRipHz muaBinnedHz ratio sptDays start stop MuaZbyDay DayT;
      
end

cat.catRip = catRip;

% TL useful to have a field that says what index of cat.dayT each ripple is
% within

end













