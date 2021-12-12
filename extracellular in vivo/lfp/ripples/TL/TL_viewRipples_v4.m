function TL_viewRipples_v4(ripples , expInfo , normalize , cr , sortVar , sortDir , bursts)

%% Tool to view data from a single recording (not catMouse output)
% Runs on output of TL_ripDetect

% INPUTS ->

% [data] : ripple structure output from TL_ripDetect
% [n] : 1 to normalize axes across all data, 0 to not
% [currentRip] : index of ripple to display initially
% [sortVar] : variable to sort the data by.
% options : 'time' , 'length' , 'power' , 'spikes'
% [sortdir] : 'ascend' or 'descend'
% [bursts] : look at bursts of data instead of single data

%% global vars
% currentRip = [];
% data = [];
% indx = [];

% muH = [];
clear n data currentRip ax aspRat scatterData;

global n data ax currentRip aspRat axSpkt ylimz indx muH scatterData artList;

currentRip = cr; clear cr;
ax= [];
ylimz = [];
data = [];
indx = [];
muH = [];
artList = [];
n = normalize;
scatterData = [];
data = ripples.detect;

%% Bin mu spiketimes for histogram (5 ms bins)
for d = 1 : length(data.spkT)
    
    %closest bins to 5 ms (assuming now 1250k samprate for lfp)
    rL = data.ripLog{d};
    bins = [0:6:length(rL)]/1250 - (find(rL,1) - 1)/1250;
    m = rem(length(rL),6);
    if m ~= 0
        bins(end+1) = (length(rL))/1250 - (find(rL,1) - 1)/1250;
    end
    h = histcounts(data.spkT{d} , bins);
    hr = repelem(h(1:end-1),6);
    if m ~= 0
        hr = [hr , repelem(h(end),m)];
    else
        hr = [hr , repelem(h(end),6)];
    end
    
    muH{d,1} = hr;
    clear h hr RL bins m;
    
end
%% Arrange data

fn = fieldnames(data);

if strcmp(sortVar , 'time')
    [~ , indx] = sort(data.lfpIndx(: , 1) , sortDir);
    for f = 1 : length(fn)
        data.(fn{f}) = data.(fn{f})(indx,:);
    end
    muH = [muH(indx)];
end

if strcmp(sortVar , 'length')
    dur = data.lfpIndx(:,2) - data.lfpIndx(:,1);
    [~ , indx] = sort(dur , sortDir);
    for f = 1 : length(fn)
        data.(fn{f}) = data.(fn{f})(indx,:);
    end
    muH = [muH(indx)];
end

if strcmp(sortVar , 'power')
    %     r = cell2mat(cellfun(@(z) max(z.^2) , ripples.detect.filtered , 'uniformoutput' , false));
    for p = 1 : size(data.lfpIndx,1)
        pw(p) = max(data.envelope{p}(logical(data.ripLog{p})));
    end
    [~ , indx] = sort(pw , sortDir);
    for f = 1 : length(fn)
        data.(fn{f}) = data.(fn{f})(indx,:);
    end
    muH = [muH(indx)];
end

if strcmp(sortVar , 'spikes')
    % calculate evoked hz
    bMs = 1000 * (find(data.ripLog{1} , 1) - 1) / 1250;
    ripS = (data.lfpIndx(:,2) - data.lfpIndx(:,1)) / 1250;
    bSpkHz = cell2mat(cellfun(@(z) sum(z < 0) , data.spkT , 'uniformoutput' , false)) / (bMs/1000) ;
    for r = 1 : length(ripS)
        evkHz(r) = sum(data.spkT{r} >= 0 & data.spkT{r} < ripS(r)) / ripS(r);
    end
    [~ , indx] = sort(evkHz , sortDir);
    for f = 1 : length(fn)
        data.(fn{f}) = data.(fn{f})(indx,:);
    end
    muH = [muH(indx)];
end

%% set up figure and axes

% 1 inch per 30 ms ripple
aspRat = 1250*30/1000;

fig = figure('units' , 'inches' , 'position' , [1 0.5 8 8] , ...
    'visible' , 'on' , 'name' , [ripples.expInfo.subject{1} ' ' num2str(ripples.expInfo.date{1})] , ...
    'KeyPressFcn' , @keyPress);

% axes for unfiltered trace
ax.unfiltered = axes('parent', fig , 'units' , 'inches' , 'position' , [4.5 6.5 diff(data.lfpIndx(currentRip,:))/aspRat 1] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 10);
hold on;
ax.unfiltered.YLabel.String = 'Raw Trace (uV)';

% Make scale bar..
scaleb = axes('parent' , fig , 'units' , 'inches' , 'position' , [4.75 6.4 15/aspRat 1] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'none' , 'XColor' , 'k' , 'XLim' , [0 15] , 'XTick' , [0:5:15] , 'XTicklabel' , {} , 'fontname' , 'arial' , 'fontsize' , 10);
scaleb.XLabel.String = '15 ms';

% axes for ripple filtered trace
ax.filtered = axes('parent', fig , 'units' , 'inches' , 'position' , [4.5 5 diff(data.lfpIndx(currentRip,:))/aspRat 1] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 10);
hold on;
ax.filtered.YLabel.String = 'Ripple Band (z)';

% axes for envelope
% ax.envelope = axes('parent', fig , 'units' , 'inches' , 'position' , [1 3.5 diff(data.lfpIndx(1,:))/aspRat 1] , 'box' , 'off' , ...
%     'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 10);
% hold on;
%
% % axes for zMua
ax.zMua = axes('parent', fig , 'units' , 'inches' , 'position' , [4.5 3.5 diff(data.lfpIndx(currentRip,:))/aspRat 1] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 10);
hold on;
ax.zMua.YLabel.String = 'Multi-Unit Firing (z)';

% axes for multiunit spike histogram
axSpkt = axes('parent', fig , 'units' , 'inches' , 'position' , [4.5 2 diff(data.lfpIndx(currentRip,:))/aspRat 1] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'k' , 'XColor' , 'none' , 'fontname' , 'arial' , 'fontsize' , 10);
hold on;
axSpkt.YLabel.String = '# MU Spikes';

artButt = uicontrol('parent' , fig , 'style' , 'pushbutton' , ...
    'units' , 'inches' , 'position' , [1 6.5 , 2.5 0.25] , ...
    'fontname' , 'arial' , 'fontsize' , 10 , 'String' , 'Remove Artifacts and Close' , ...
    'callback' , {@cb_artButt , ripples , expInfo});





% Scatter plot axes

scatAx = axes('parent' , fig , 'units' , 'inches' , 'position' , [0.5 , 2.75 , 3 , 3] , 'box' , 'off' , ...
    'color' , 'none' , 'YColor' , 'k' , 'fontname' , 'arial' , 'fontsize' , 10)


%% Make scatter data
% to get true time, need to see if sample is in a separate block, and find
% the start time of that block and add it.
ns = cell2mat(expInfo.nsec);
st = cell2mat(cellfun(@(z) datenum(z) , expInfo.startTime , 'uniformoutput' , false));
st = st - st(1);
allT = [];
for s = 1 : length(ns)
    i = data.lfpIndx(:,1)/1250 < sum(ns(1:s));
    currT = data.lfpIndx(i,1)/1250;
    currT = currT + st(s);
    allT = [allT ; currT];
    clear currT;
end

% find mean and peak zMua for each ripple
for a = 1 : length(allT)
    dz = data.zMua{a}(logical(data.ripLog{a}));
    ez = data.envelope{a}(logical(data.ripLog{a}));
    mnZ(a) = mean(dz);
    mxZ(a) = max(dz);
    mnEZ(a) = mean(ez);
    mxEZ(a) = max(ez);
end

% TL best will be to have the ripple file have the true times locked to t =
% 0 and the block that each ripple is found in (to find the gap between
% blocks will need expInfo;







% find ylims for each type of plot
if n
    fn = fieldnames(ax)
    for f = 1 : length(fn)
        d = data.(fn{f});
        ma = max(cell2mat(cellfun(@(z) max(z) , d , 'uniformoutput' , false)));
        mi = min(cell2mat(cellfun(@(z) min(z) , d , 'uniformoutput' , false)));
        di = ma - mi;
        ylimz.(fn{f}) = [mi - 0.05 * di , ma + 0.05 * di];
    end
end
makePlots;

    function cb_artButt(~ , ~ , ripples , expInfo)
        global artList
        fn = fieldnames(ripples.detect)
        for f = 1 : length(fn)
            ripples.detect.(fn{f})(artList,:) = []
        end
        
        i = find(cell2mat(expInfo.set) == ripples.set , 1);
        f = strfind(expInfo.basepath{i} , 'Dates') - 1
        savePath = [expInfo.basepath{i}(1:f) , 'Subjects' , filesep , expInfo.subject{i} , filesep , num2str(expInfo.date{i}) , filesep , 'ripples_' num2str(ripples.set) '.mat'];
save(savePath , 'ripples');
close



    function currentRip = keyPress(a , b)
        global ax data indx aspRat currentRip axSpkt artList
        
        if strcmp(b.Key , 'leftarrow') | strcmp(b.Key , 'rightarrow')
            if strcmp(b.Key , 'leftarrow') & currentRip > 1
                currentRip = currentRip - 1;
            end
            
            if strcmp(b.Key , 'rightarrow') & currentRip < size(data.lfpIndx , 1)
                currentRip = currentRip + 1;
            end
            
            fn = fieldnames(ax);
            for f = 1 : length(fn)
                c = get(ax.(fn{f}) , 'children');
                for cc = 1 : length(c)
                    delete(c(cc));
                end
            end
            
            c = get(axSpkt , 'children');
            for cc = 1 : length(c)
                delete(c(cc));
            end
        else
            if strcmp(b.Key , 'a')
                artList = [artList; indx(currentRip)];
            else if strcmp(b.Key , 'r')
                    
                    artList = artList(~ismember(artList , indx(currentRip)));
                end
            end
        end
    

makePlots;%(ax , data , currentRip , aspRat);

    function makePlots
        global ax data aspRat currentRip n ylimz indx muH axSpkt;
        fn = fieldnames(ax);
        for f = 1 : length(fn)
            d = data.(fn{f}){currentRip};
            L = logical(data.ripLog{currentRip});
            %     if strcmp(fn{f} , 'unfiltered')
            %         L = repelem(L , 3 , 1);
            %     end
            xval = 1 : length(d);
            plot(ax.(fn{f}) , xval , d , 'linestyle' , '-' , 'color' , [0.4 0.4 0.4] , 'linewidth' , 0.25 , 'marker' , 'none');
            if ~strcmp(fn{f} , 'unfiltered')
                plot(ax.(fn{f}) , xval(L) , d(L) , 'linestyle' , '-' , 'color' , 'k' , 'linewidth' , 0.5 , 'marker' , 'none');
            end
            ax.(fn{f}).XLim = [0 , length(d)];
            ax.(fn{f}).Position(3) = diff(data.lfpIndx(currentRip,:))/aspRat;
            
            if strcmp(fn{f} , 'filtered')% | strcmp(fn{f} , 'unfiltered')
                d = data.envelope{currentRip};
                
                p = plot(ax.(fn{f}) , xval(L) , d(L) , 'linestyle' , '-' , 'color' , [1 0.1 0.1]  , 'linewidth' , 1 , 'marker' , 'none');
                uistack(p,'bottom');
                p = plot(ax.(fn{f}) , xval , d , 'linestyle' , '-' , 'color' , [1 0.1 0.1] , 'linewidth' , 0.25 , 'marker' , 'none');
                uistack(p,'bottom');
            end
            if n
                ax.(fn{f}).YLim = ylimz.(fn{f});
            end
            
            if f == 1
                ax.(fn{f}).Title.String = [num2str(currentRip) ' ' num2str(indx(currentRip))];
            end
        end
        
        bar(axSpkt , xval , muH{currentRip} , 1 , 'facecolor' , 'k' , 'linestyle' , 'none');
        plot(axSpkt , [repmat(find(diff(data.ripLog{currentRip}) == 1)+1,1,2)] , axSpkt.YLim , 'color' , 'r' , 'linestyle' , ':' , 'linewidth' , 0.75);
        plot(axSpkt , [repmat(find(diff(data.ripLog{currentRip}) == -1),1,2)] , axSpkt.YLim , 'color' , 'r' , 'linestyle' , ':' , 'linewidth' , 0.75);
        c = sum(data.spkT{currentRip} >= 0 & data.spkT{currentRip} < sum(data.ripLog{currentRip})/1250);
        text(axSpkt , find(diff(data.ripLog{currentRip}) == -1) , axSpkt.YLim(2) , num2str(c) , 'color' , 'b' , 'fontname' , 'arial' , 'fontsize' , 10 , 'horizontalalignment' , 'right');
        axSpkt.Position(3) = diff(data.lfpIndx(currentRip,:))/aspRat;
        
        
        function cb_artifact
            
            global indx currentRip artList
            
            
            %zMua plot
            % ax.mua.zMua =
            
            % fn = fieldnames(ax.mua)
            % for f = 1 : length(fn)
            % ax.mua.(fn{f}).Position(3) = diff(data.lfpIndx(currentRip,:))/aspRat;
            % end
