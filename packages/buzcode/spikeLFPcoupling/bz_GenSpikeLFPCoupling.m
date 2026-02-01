function [SpikeLFPCoupling] = bz_GenSpikeLFPCoupling(spikes,LFP,varargin)
% SpikeLFPCoupling = GenSpikeLFPCoupling(spikes,LFP)
%
%INPUT
%   spikes          structure with fields (from bz_getSpikes)
%                           .times
%                       -or-
%                       {Ncells} cell array of spiketimes for each cell
%   LFP                 structure with fields (from bz_GetLFP)
%                           lfp.data
%                           lfp.timestamps
%                           lfp.samplingRate
%   (optional)
%       'frange'
%       'tbin'
%       'waveparms'...
%       'nfreqs'
%       'int'
%       'ISIpower'      true: calculate mutual information between power
%                       ISI distribution
%       'cellLFPchannel'  The local LFP channel associated with each cell (if
%                       sites are in different electroanatomical regions)
%       'cellsort'        'pca','none','sortf','sortorder','celltype','rate'
%       'controls'        'thinspikes','jitterspikes','shufflespikes'
%       'downsample'
%       'cellclass'     cell array with string label for class of each cell
%       'jittersig'     true/false for jittered spikes significance
%       'showFig'       true/false
%       'spikeLim'      limit number of spikes to look at for each cell
%                       (randomly omits spikes, default: Inf)
%
%OUTPUT
%        SpikeLFPCoupling
%           .freqs
%           .pop
%               .popcellind
%               .cellpopidx
%           .cell                all results are [cells x freqs x channels]
%               .ratepowercorr
%               .ratepowersig
%               .spikephasemag
%               .spikephaseangle
%               .spikephasesig
%
%           .detectorinfo
%               .detectorname
%               .detectiondate
%               .detectionintervals
%               .detectionchanne;
%
%TO DO
%   -Put synch/spikerate in same section
%   -multiple ints - this (along with controls) will require making
%   subfunctions to loop
%   -clean and buzcode
%
%DLevenstein 2016

%% inputParse for Optional Inputs and Defaults
p = inputParser;


checkInt = @(x) size(x,2)==2 && isnumeric(x) || isa(x,'intervalSet');
checkFrange = @(x) isnumeric(x) && length(x(1,:)) == 2 && length(x(:,1)) == 1;

validSorttypes = {'pca','none','sortf','sortorder','celltype','rate'};
checkSorttype = @(x) any(validatestring(x,validSorttypes)) || size(x) == [1,2];

addParameter(p,'int',[0 Inf],checkInt)
addParameter(p,'frange',[1 128],checkFrange)
addParameter(p,'nfreqs',100,@isnumeric)
addParameter(p,'ncyc',7,@isnumeric)
addParameter(p,'synchdt',0.005,@isnumeric)
addParameter(p,'synchwin',0.02,@isnumeric)
addParameter(p,'sorttype','rate',checkSorttype)
addParameter(p,'DOWNSAMPLE',false,@isnumeric)
addParameter(p,'cellclass',0)
addParameter(p,'channel',[])
addParameter(p,'jittersig',false)
addParameter(p,'showFig',true)
addParameter(p,'saveFig',false)
addParameter(p,'saveMat',false)
addParameter(p,'ISIpower',false)
addParameter(p,'spikeLim',Inf)

parse(p,varargin{:})
%Clean up this junk...
int = p.Results.int;
nfreqs = p.Results.nfreqs;
frange = p.Results.frange;
ncyc = p.Results.ncyc;
synchdt = p.Results.synchdt;
jittersig = p.Results.jittersig;
SHOWFIG = p.Results.showFig;
SAVEMAT = p.Results.saveMat;
figfolder = p.Results.saveFig;
usechannel = p.Results.channel;
ISIpower = p.Results.ISIpower;
spikeLim = p.Results.spikeLim;
subpop = p.Results.cellclass;

%% Deal with input types

if isa(spikes,'cell')
    spikes_temp.numcells = length(spikes);
    spikes_temp.times = spikes;
    spikes = spikes_temp;
    clear spikes_temp
end

if isa(int,'intervalSet')
    int = [Start(int,'s'), End(int,'s')];
end


if ~isempty(usechannel)
    usechannel = ismember(LFP.channels,usechannel);
    LFP.data = LFP.data(:,usechannel);
end

%Downsampling
if p.Results.DOWNSAMPLE
    assert((LFP.samplingRate/p.Results.DOWNSAMPLE)>2*max(frange),'downsample factor is too big...')
    [ LFP ] = bz_DownsampleLFP(LFP,p.Results.DOWNSAMPLE);
end
%% Subpopulations
if isequal(subpop,0)
    numpop = 1;
    popcellind = {1:length(spikes.times)};
    pops = {'rs'};
elseif isequal(subpop,'done')
else
    emptycells =cellfun(@isempty,subpop);
    pops = unique(subpop(~emptycells));
    subpop(emptycells) = {'nopop'};
    numpop = length(pops);
    for pp = 1:numpop
        popcellind{pp} = find(ismember(subpop,pops{pp}));
    end
end

cellpopidx = zeros(1,spikes.numcells);


%% Calculate spike matrix
spikemat = bz_SpktToSpkmat(spikes,'binsize',p.Results.synchwin,'dt',synchdt);

inint = InIntervals(spikemat.timestamps,int);
spikemat.data = spikemat.data(inint,:);
spikemat.timestamps = spikemat.timestamps(inint);


%% Calculate ISIs (if needed), Apply interval/spike number limits here

%Take only spike times in intervals
spikes.inint = cellfun(@(X) InIntervals(X,int),spikes.times,'UniformOutput',false);

% if ISIpower
%     %Calculate (pre/postceding ISI) for each spike
%     spikes.ISIs_prev = cellfun(@(X) [nan; diff(X)],spikes.times,'UniformOutput',false);
%     spikes.ISIs_next = cellfun(@(X) [diff(X); nan],spikes.times,'UniformOutput', false);
% end

%Apply spikelimit here
spikes.toomany = cellfun(@(X) (sum(X)-spikeLim).*((sum(X)-spikeLim)>0),spikes.inint);
spikes.numcells = length(spikes.times);
for cc = 1:spikes.numcells
    spikes.inint{cc}(randsample(find(spikes.inint{cc}),spikes.toomany(cc)))=false;
end

%Restrict to up to spikelim spikes in interval
spikes.times = cellfun(@(X,Y) X(Y),spikes.times,spikes.inint,'UniformOutput',false);

%% Processing LFP - filter etc

%HERE: loop channels
for cc = 1:length(LFP.channels)
    chanID = LFP.channels(cc);
    if length(LFP.channels)>1
        fprintf('Channel %d (%d of %d)\n',chanID,cc,length(LFP.channels))
    end
    
    switch nfreqs
        case 1
            %Single frequency band - filter/hilbert
            LFP_filt = bz_Filter(LFP,'passband',frange,'order',ncyc,...
                'intervals',int,'filter','fir1','channels',chanID,'fast',false);
            LFP_filt.data = LFP_filt.hilb;
            freqs = [];
            clear filtered
            
            %Normalize Power to Mean Power
            LFP_filt.data = LFP_filt.data./mean(abs(LFP_filt.data));
            
        otherwise
            %Multiple frequencies: Wavelet Transform
            LFP_filt = bz_WaveSpec(LFP,'intervals',int,'showprogress',true,...
                'ncyc',ncyc,'nfreqs',nfreqs,'frange',frange,'chanID',chanID);
            freqs = LFP_filt.freqs;
            
            %Normalize power to mean power for each frequency
            LFP_filt.data = bsxfun(@(X,Y) X./Y,LFP_filt.data,nanmean(abs(LFP_filt.data),1));
    end
    
    %Get Power/Phase at each spike matrix time point and each spike time
    % spike mat is a spike count vector for each cell calculated in dt
    % steps with binsize overlap. the overlap is done with a moving
    % average. here, the filtered power in each frequency band is retrieved
    % for each spike count bin.
    spikemat.filtLFP = interp1(LFP_filt.timestamps,LFP_filt.data,spikemat.timestamps,'nearest');
    
    
    %Get complex-valued filtered LFP at each spike time
    % here, the filtered power in each frequency band is retrieved
    % for each spike of each unit
    if spikes.numcells>50
        disp('Interpolating LFP at each spike... If this is prohibitive (time or RAM), try using ''spikeLim''')
    end
    for nn = 1:spikes.numcells
        bz_Counter(nn,spikes.numcells,'Interpolating Cell')
        spikes.filtLFP{nn} = interp1(LFP_filt.timestamps,LFP_filt.data,spikes.times{nn},'nearest');
    end
    
    
    %% Cell Rate-Power Modulation
    
    %Spike-Power Coupling
    % here you get the correlation between the number of spikes in each dt
    % bin and the power of the lfp in each frequency band. the output is a
    % 2d matrix. it seems to me that the step size of spike counts should
    % be adjusted according to the frequency of interest. note that the filtered data is the hilbert
    % transform and here the amplitude is taken
    [ratepowercorr(:,:,cc),ratepowersig(:,:,cc)] = corr(spikemat.data,abs(spikemat.filtLFP),'type','spearman','rows','complete');
    
    %% Cell Spike-Phase Coupling
    
    %Find filtered LFP at the closest LFP timepoint to each spike.
    
    nunits = length(spikes.times);
    numcells = nunits;
    for iunit = 1 : nunits
        spkZ = mean(spikes.filtLFP{iunit}, 1, 'omitnan');
        spkMag(iunit, :) = abs(spkZ);
        spkPhase(iunit, :) = angle(spkZ);
    end
    spikephasemag = spkMag;
    spikephaseangle = spkPhase;
    
end
% Example Figure : Phase-Coupling
% if phmag >0.1
%     figure
%             rose(angle(spkLFP_fn(:, 3)))
%            % set(gca,'ytick',[])
%             hold on
%              polar(angle(spkLFP_fn(:, 3)),abs(spkLFP_fn(:, 3)).*40,'k.')
%             % hold on
%              %compass([0 phangle],[0 phmag],'r')
%              compass(rvect.*700,'r')
%         delete(findall(gcf,'type','text'));
%         % delete the text objects

% this takes an insane amount of time
if jittersig
    %Jitter for Significane
    numjitt = 1;
    jitterwin = 2/frange(1);
    
    jitterbuffer = zeros(numcells,nfreqs,numjitt);
    for jj = 1:numjitt
        if mod(jj,10) == 1
            display(['Jitter ',num2str(jj),' of ',num2str(numjitt)])
        end
        jitterspikes = bz_JitterSpiketimes(spikes.times,jitterwin);
        jitterLFP = cellfun(@(X) interp1(LFP_filt.timestamps,LFP_filt.data,X,'nearest'),...
            jitterspikes,'UniformOutput',false);
        
        for iunit = 1 : nunits
            spkZ = mean(jitterLFP{iunit}, 1, 'omitnan');
            jitterbuffer(iunit, :, jj) = abs(spkZ);
        end
        
    end
    jittermean = mean(jitterbuffer,3);
    jitterstd = std(jitterbuffer,[],3);
    spikephasesig(:,:,cc) = (spikephasemag-jittermean)./jitterstd;
end


% Population Synchrony: Phase Coupling and Rate Modulation
% here is where the separation of rs and fs is useful
for pp = 1:numpop
    if length(popcellind{pp}) == 0
        popcoupling.(pops{pp}).powercorr = [];
        popcoupling.(pops{pp}).phasemag = [];
        popcoupling.(pops{pp}).phaseangle = [];
        numpop = numpop-1;
        continue
    end
    cellpopidx(popcellind{pp}) = pp;
    numpopcells = length(popcellind{pp});
    % here the mean activity across the population is calculated. not
    % sure why convert the spike counts to binary instead of averaging
    % the counts? basically gets the average number of cells that were
    % active in each time bin - a measurment of synchorny
    popsynch = sum(spikemat.data(:, popcellind{pp}) > 0, 2) ./ numpopcells;
    % standardize  to the mean, probably so that the correlation
    % considers the overall degree of synchrony
    popsynch = popsynch ./ mean(popsynch);
    
    %Calculate Synchrony-Power Coupling as correlation between synchrony and power
    [popcoupling.(pops{pp}).powercorr(:,cc)] = corr(popsynch,abs(spikemat.filtLFP),'type','spearman','rows','complete');
    
    %Synchrony-Phase Coupling (magnitude/angle of power-weighted mean resultant vector)
    %Note: need a way ehere to account for phase/power occupancy...
    %could not find a reference for what is done here
    resultvect = nanmean(abs(spikemat.filtLFP).*bsxfun(@(popmag,ang) popmag.*exp(1i.*ang),...
        popsynch,angle(spikemat.filtLFP)),1);
    
    lfpmag = abs(spikemat.filtLFP);
    lfpphase = angle(spikemat.filtLFP);
    resultvect = mean(lfpmag .* (popsynch .* exp(1i .* lfpphase)),...
        1, 'omitnan');
    
    
    
    popcoupling.(pops{pp}).phasemag(:,cc) = abs(resultvect);
    popcoupling.(pops{pp}).phaseangle(:,cc) = angle(resultvect);
    
    %         if ISIpower
    %             popcoupling.(pops{pp}).ISIpowermodulation(:,cc) = nanmean(totmutXPow(popcellind{pp},:,cc),1);
    %         end
end


clear LFP_filt

%% Output

SpikeLFPCoupling.freqs = freqs;
SpikeLFPCoupling.pop = popcoupling;
%SpikeLFPCoupling.pop.popcellind = popcellind; %update to match buzcode cell class input
%SpikeLFPCoupling.pop.cellpopidx = cellpopidx;
SpikeLFPCoupling.cell.ratepowercorr = ratepowercorr;
SpikeLFPCoupling.cell.ratepowersig = ratepowersig;
SpikeLFPCoupling.cell.spikephasemag = spikephasemag;
SpikeLFPCoupling.cell.spikephaseangle = spikephaseangle;
if ISIpower
    SpikeLFPCoupling.cell.ISIpowermodulation = totmutXPow;
end
if exist('spikephasesig','var')
    SpikeLFPCoupling.cell.spikephasesig = spikephasesig;
end

SpikeLFPCoupling.detectorinfo.detectorname = 'bz_GenSpikeLFPCoupling';
SpikeLFPCoupling.detectorinfo.detectiondate = datetime('today');
SpikeLFPCoupling.detectorinfo.detectionintervals = int;
SpikeLFPCoupling.detectorinfo.detectionchannel = LFP.channels;

if SAVEMAT
    [~, basename] = fileparts(pwd);
    savefile = fullfile(pwd, [basename, '.spkLFPcoupling.mat']);
    save(savefile,'SpikeLFPCoupling')
end

%% Figures
if SHOWFIG
    
    %Sorting (and other plot-related things)
    switch p.Results.sorttype
        case 'pca'
            [~,SCORE,~,~,EXP] = pca(ratepowercorr);
            [~,spikepowersort] = sort(SCORE(:,1));
            
            [~,SCORE,~,~,EXP_phase] = pca(spikephasemag);
            [~,spikephasesort] = sort(SCORE(:,1));
            sortname = 'PC1';
        case 'none'
            spikepowersort = 1:numcells;
        case 'fsort'
            fidx = interp1(freqs,1:nfreqs,sortf,'nearest');
            [~,spikepowersort] = sort(ratepowercorr(:,fidx));
            [~,spikephasesort] = sort(spikephasemag(:,fidx));
            sortname = [num2str(sortf) 'Hz Magnitude'];
        case 'rate'
            spkrt = cellfun(@length,spikes.times);
            [~,spikepowersort] = sort(spkrt);
            spikephasesort = spikepowersort;
            sortname = 'Firing Rate';
        otherwise
            spikepowersort = sortidx;
    end
    
    
    if length(LFP.channels)==1
        
        switch nfreqs
            case 1
                %% Figure: 1 Freq Band
                figure
                subplot(3,2,1)
                hold on
                hist(abs(ratepowercorr))
                for pp = 1:numpop
                    plot([popcoupling.(pops{pp}).powercorr,popcoupling.(pops{pp}).powercorr],...
                        get(gca,'ylim'))
                end
                xlabel('Rate-Power Correlation')
                title('Rate-Power Coupling')
                ylabel('# Cells')
                subplot(3,2,[2 4])
                for pp = 1:numpop
                    polar(spikephaseangle(popcellind{pp}),(spikephasemag(popcellind{pp})),'o')
                    hold on
                    polar([0 popcoupling.(pops{pp}).phaseangle],[0 popcoupling.(pops{pp}).phasemag])
                    
                end
                title('Spike-Phase Coupling')
                subplot(3,2,3)
                for pp = 1:numpop
                    hold on
                    plot(ratepowercorr(popcellind{pp}),(spikephasemag(popcellind{pp})),'o')
                    plot(popcoupling.(pops{pp}).powercorr,popcoupling.(pops{pp}).phasemag,'*')
                end
                %LogScale('y',10)
                xlabel('Rate-Power Correlation')
                ylabel('Spike-Phase Coupling Magnitude')
                
                if figfolder
                    NiceSave('SpikeLFPCoupling',figfolder,[])
                end
                
            otherwise
                %% Figure: Spectrum
                posnegcolor = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
                
                figure
                subplot(3,2,1)
                hold on
                for pp = 1:numpop
                    plot(log2(freqs),popcoupling.(pops{pp}).powercorr)
                end
                plot(log2(freqs([1 end])),[0 0],'k--')
                axis tight
                box off
                LogScale('x',2)
                title('Pop. Synchrony - Power Correlation')
                xlabel('f (Hz)');ylabel('rho')
                
                %     subplot(3,2,3)
                %         hold on
                %             [hAx,hLine1,hLine2] = plotyy(log2(freqs),cat(1,synchcoupling.phasemag),...
                %                 log2(freqs),mod(cat(1,synchcoupling.phaseangle),2*pi))
                %         LogScale('x',2)
                %         hLine1.LineStyle = 'none';
                %         hLine2.LineStyle = 'none';
                %         hLine1.Marker = 'o';
                %         hLine2.Marker = '.';
                %         title('Pop. Synchrony - Phase Coupling')
                %         xlabel('f (Hz)');
                %         ylabel(hAx(1),'Phase Coupling Magnitude')
                %         ylabel(hAx(2),'Phase Coupling Angle')
                
                subplot(3,2,3)
                hold on
                for pp = 1:numpop
                    plot(log2(freqs),popcoupling.(pops{pp}).phasemag)
                end
                axis tight
                box off
                LogScale('x',2)
                title('Pop. Synchrony - Phase Coupling')
                xlabel('f (Hz)');
                ylabel('Phase Coupling Magnitude')
                subplot(3,2,5)
                hold on
                for pp = 1:numpop
                    plot(log2(freqs),popcoupling.(pops{pp}).phaseangle,'o')
                    plot(log2(freqs),popcoupling.(pops{pp}).phaseangle+2*pi,'o')
                end
                LogScale('x',2)
                title('Pop. Synchrony - Phase Coupling')
                xlabel('f (Hz)');
                axis tight
                box off
                ylim([-pi 3*pi])
                ylabel('Phase Coupling Angle')
                
                
                subplot(2,2,2)
                imagesc(log2(freqs),1:spikes.numcells,real(ratepowercorr(spikepowersort,:)))
                xlabel('f (Hz)');ylabel(['Cell - Sorted by ',sortname])
                LogScale('x',2)
                title('Rate - Power Correlation')
                colormap(gca,posnegcolor)
                ColorbarWithAxis([-0.2 0.2],'rho')
                subplot(2,2,4)
                imagesc(log2(freqs),1:spikes.numcells,spikephasemag(spikephasesort,:))
                xlabel('f (Hz)');ylabel(['Cell - Sorted by ',sortname])
                title('Spike - Phase Coupling')
                LogScale('x',2)
                caxis([0 0.2])
                colorbar
                
                if figfolder
                    NiceSave('SpikeLFPCoupling',figfolder,[])
                end
                
        end
        
        
        
    end
end

end







%     %% Spike-Phase Coupling function
%     %takes spike times from a single cell and caluclates phase coupling magnitude/angle
%     function [phmag,phangle] = spkphase(spkLFP_fn)
%         %Spike Times have to be column vector
%             if isrow(spkLFP_fn); spkLFP_fn=spkLFP_fn'; end
%             if isempty(spkLFP_fn); phmag=nan;phangle=nan; return; end
%         %Calculate (power normalized) resultant vector
%         rvect = nanmean(abs(spkLFP_fn).*exp(1i.*angle(spkLFP_fn)),1);
%         phmag = abs(rvect);
%         phangle = angle(rvect);
%
%         %% Example Figure : Phase-Coupling
%         % if phmag >0.1
%         % figure
%         %         rose(phase4spike)
%         %        % set(gca,'ytick',[])
%         %         hold on
%         %          polar(phase4spike,power4spikephase.*40,'k.')
%         %         % hold on
%         %          %compass([0 phangle],[0 phmag],'r')
%         %          compass(rvect.*700,'r')
%         %     delete(findall(gcf,'type','text'));
%         %     % delete the text objects
%         % end
%     end



