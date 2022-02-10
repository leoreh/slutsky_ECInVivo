
% idetifying same cells from different contineous recordings

% STEPS
% (1) calucalte isolation distance between all clusters from two
% adjacent recordings and combine similar clusters. 
% (2) spike detection on the last / first x hours of two adjacent
% recordings. NO NEED FOR SPIKE DETECTION. just take all spikes from all
% clusters, ignore their identiy and project them on feature space. then
% move to step 3. 
% (3) k-means hard clustering of the spikes from step (2) using the
% cluster centers obtained in step (1)
% (4) examine CCs between the clusters in step (3) and all the clusters
% from the original recordings. two clusters from the original recordings
% that have a (very) significant peak at 0 ms time lag with one of the new
% clusters belong to the same unit
% (5) validate waveform similarity and ACGs manually 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

margintime = 1 * 60 * 60;   % time to use from edges of recordings [s]
sgrps = [1 : 1];            % selected groups to include in the analysis 
npca = 3;                   % number of PCs per spk grp
graphics = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1 - load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path to recording sessions. must be ordered chornically
basepaths = {'D:\Google Drive\Data\lh96\lh96_211201_070100',...
    'D:\Google Drive\Data\lh96\lh96_211202_070500'};

% initialize
lastclu = 0;
fet = cell(1, 2);
clu = cell(1, 2);
spktimes = cell(1, 2);
for ifile = 1 : 2
    
    % session params
    session = CE_sessionTemplate(basepaths{ifile}, 'viaGUI', false,...
        'forceDef', true, 'forceL', false, 'saveVar', true);
    spkgrp = session.extracellular.spikeGroups.channels;
    fs = session.extracellular.sr;
    marginsamples = margintime * fs;
    
    for igrp = sgrps
        
        % load data
        grpch = spkgrp{igrp};
        nfet = length(grpch) * npca;       
        fetTmp = loadNS('datatype', 'fet', 'session', session, 'grpid', igrp,...
            'nfet', nfet + length(grpch) + 1);
        cluTmp = loadNS('datatype', 'clu', 'session', session, 'grpid', igrp);
        
        % clip spikes to last margintime
        if ifile == 1
            lastspk = max(fetTmp(end, :));
            spkidx = fetTmp(:, end) > lastspk - marginsamples & fetTmp(:, end) < Inf;
        elseif ifile == 2
            spkidx = fetTmp(:, end) > 0 & fetTmp(:, end) < marginsamples;
        end                      
        
        % remove noise clusters
        spkidx = spkidx & cluTmp ~= 0 & cluTmp ~= 1;
        
        % arrange in cell
        fet{ifile} = [fet{ifile}; fetTmp(spkidx, 1 : nfet)];
        clu{ifile} = [clu{ifile}; cluTmp(spkidx) + lastclu];
        spktimes{ifile} = [spktimes{ifile}; fetTmp(spkidx, end)];

        lastclu = max(unique(cluTmp));

        % clear memory
        clear fetTmp spkTmp cluTmp
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2 - merge similar clusters from the two recordings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UPDATE 05.02.22 - realized logical mistake - spikes include both the
% timing and waveoform so it is impossible to use one as a ground truth for
% the other, because naturally if clusters are megred based on waveform
% similairty, their timings will be the same too. so there is no way to
% avoid fancy statistics. but maybe if re-clustering is done by waveform in
% feature space, and the test is done on CCs / ACs w/ bootstraping, the
% validation will be strong enough. i.e., recluster the spikes from
% recording A to the centroids of recording B and examine the similarity in
% CCs. then repeat a 1000 times but each to recluster the spikes to random
% centroids. perhaps this is still just a fancy way to use spike timing as
% a validation to waveform similairty. perhaps there is a way to strip the
% original spike from their timing, replace it with simulated timings
% according to the clusters ac, and then see if the ac similiarity between
% the new and original cluster are similar. I think i've got - recluster A
% (only) according to centroids of B, and test the similarity between ACs
% of B and the new clusters. Because only the spikes of A are in the new
% clusters, AC similarity contains totally new information (i.e. there is
% no pre-determined bias that similar waveforms will elicit similiar ACs).
% it is still not bullet proof to FP (e.g. from different mice), but its
% statistical power include both the chance that two waveforms are similar
% AND the chance that their firing patterns are similar. if the
% re-clustering and constructino of ACs is done on a fixed amount of time
% rather then a fixed amount of spikes, than the non-normalized ACs (in
% counts) will also include firing rate alongside pattern. before the ACs
% step, the re-clustering must be validated as successful only for clusters
% that maintained x% of their original spikes. 

% DIFFERENT APPROACH. re-assign spikes from the second recording to
% cluster centers from the first recording. then for each new cluster
% check if the "reshuffling by waveform" maintained a significnat 0 lag
% peak and if so, if the new and original waveforms are similar.  

% if clusters are not merged based on waveform, the "new" spikes will be
% arbiturary assigned to each one of them and thus may greatly reduce the
% chance to detect significant correlations. this is more important then
% mistankly assigning "noise" spikes to a cluster. note here i assume
% clusters from the same session are properly isolated.

% initialize.
uclu1 = unique(clu{1});
uclu2 = unique(clu{2});
iDist = nan(length(uclu1),  length(uclu2));

% loop through pairs and calc distance. perhaps this should be done
% recursively
for u1 = 1 : length(uclu1)
    unit1idx = clu{1} == uclu1(u1);
    if sum(unit1idx) < nfet
        continue
    end
    for u2 = 1 : length(uclu2)
        unit2idx = clu{2} == uclu2(u2);
        if sum(unit2idx) < nfet
            continue
        end
        
        d2 = mahal([fet{1}(unit1idx, :); fet{2}(unit2idx, :)],...
            fet{2}(unit2idx, :));
        iDist(u1, u2) = median(d2(1 : sum(unit1idx)));
    end
end

% cell 2 mat
clu = cat(1, clu{:});
fet = cat(1, fet{:});
spktimes = cat(1, spktimes{:});

% combine similar clusters
% distThr = 20; 
% clear upair
% [upair(:, 1), upair(:, 2)] = find(iDist < distThr);
% [unpair, ia, ib] = unique(upair(:, 1));
% for ipair = 1 : length(unpair)
%     clu(clu == unpair(ipair)) = uclu2(ia(ipair));
% end


uclu = unique(clu);
nclu = length(uclu);
sclu = randperm(length(uclu), length(uclu));

% graphics ----------------------------------------------------------------
if graphics
    nplots = 6;
    maxspks = 3000;
    pcpair = randperm(nfet, nplots * 2);
    pcpair = reshape(pcpair, nplots, 2);
    
    % initialize
    fh = figure;
    cmap = colormap(jet(length(uclu)));
    cmap = cmap(randperm(length(cmap)), :);
    
    for u1 = 1 : length(sclu)
        cluidx = find(clu == sclu(u1));
        spkidx = 1 : min([maxspks, length(cluidx)]);
        for iplot = 1 : nplots
            subplot(nplots / (nplots / 2), nplots / 2, iplot)
            plot(fet(cluidx(spkidx), pcpair(iplot, 1)),...
                fet(cluidx(spkidx), pcpair(iplot, 2)), '.', 'MarkerSize', 5,...
                'Color', cmap(u1, :))
            hold on
            xlabel(sprintf('PC%d', pcpair(iplot, 1)))
            ylabel(sprintf('PC%d', pcpair(iplot, 2)))
            xticks = [];
            yticks = [];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3 - re-assign the spikes to the new clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find centroid location for each cluster
cntr = zeros(length(uclu), nfet);
for u1 = 1 : length(uclu)
    cluidx = find(clu == uclu(u1));
    cntr(u1, :) = mean(fet(cluidx, :));
end

fprintf('\nre-assigning spikes to clusters...\n')
clunew = kmeans(fet, [], 'start', cntr, 'MaxIter', 200);

uclunew = unique(clunew);
nclunew = length(uclunew);

spktimesnew = cell(nclunew, 1);
for iunit = 1 : nclunew
    cluidx = clunew == uclunew(iunit);
    spktimesnew{iunit} = spktimes(cluidx);
end

spktimescell = cell(nclu, 1);
for iunit = 1 : length(uclu)
    cluidx = clu == uclu(iunit);
    spktimescell{iunit} = spktimes(cluidx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4 - find significant correlations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

monosyn = monoSyn_wrapper('spktimes', [spktimesnew; spktimescell], 'basepath', pwd,...
    'winCalc', [0, Inf], 'saveVar', false, 'graphics', true,...
    'forceA', true, 'fs', fs, 'saveFig', false,...
    'wv', [], 'wv_std', []);

% get idx comparebale between uclu and monosyn. note order of cell
% concatination in monosyn input
idx = 1 : nclunew;
idx1 =  nclunew + 1 : nclunew + length(uclu1);
idx2 =  idx1(end) + 1 : idx1(end) + length(uclu2);

idxms = [uclunew; uclu1; uclu2];

% select threesomes 
[upair(:, 1), upair(:, 2)] = find(monosyn.eSig);
upair = sortrows(upair, 1);

% sig corr between recordings and new
new2one = ismember(upair(:, 1), idx) & ismember(upair(:, 2), idx1);
new2two = ismember(upair(:, 1), idx) & ismember(upair(:, 2), idx2);
one2two = ismember(upair(:, 1), idx1) & ismember(upair(:, 2), idx2);    % sanity check only

% threesoms. must make sure there is only a single fit between 1 and 2
upair(new2one, :)
upair(new2two, :)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 5 - manually review pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load spk for plotting the waveform or swv.


