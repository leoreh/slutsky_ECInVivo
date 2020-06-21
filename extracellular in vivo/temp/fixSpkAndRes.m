% fixClustering

% since the waveforms were detrended in ks2ns, we can assume the beginning
% and end of a wave are equal to zero and thus circhshift can be used to
% recenter the minimum.

% problem: if the waveform includes two spikes, the larger one will always
% be centered even if the cluster is refers to the smaller one.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

basepath = 'D:\dataTemp\080610';
cd(basepath)

[~, basename] = fileparts(basepath);
load([basename '.session.mat'])
spkgrp = session.extracellular.spikeGroups.channels;
ngrp = length(spkgrp);  % number of tetrodes
nsamps = 32;            % waveform length
psamp = 16;             % peak samples
dt = 8;                 % dead time [samples]. 8 @ 20000 Hz = 400 us

% spike groups to work on
if isempty(grp)
    grp = 1 : ngrp;
end

% files
clufiles = dir([basepath filesep '*.clu*']);
clunames = {clufiles.name};
resfiles = dir([basepath filesep '*.res*']);
resnames = {resfiles.name};
spkfiles = dir([basepath filesep '*.spk*']);
spknames = {spkfiles.name};
if length(resfiles) ~= length(clufiles) ||...
        length(clufiles) ~= length(spkfiles)
    error('different number of res / clu / spk files')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% go over groups
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = grp
    
    nchans = length(spkgrp{i});
    
    % load data
    fprintf('Loading data of group %d / %d\n', i, ngrp)
      
    cluname = fullfile(clufiles(i).folder, clufiles(i).name);
    fid = fopen(cluname, 'r');
    nclu = fscanf(fid, '%d', 1);
    clu = fscanf(fid, '%f', Inf);
    fclose(fid);
    cluGrp = unique(clu);
    nspks = length(clu);
    
    resname = fullfile(resfiles(i).folder, resfiles(i).name);
    fid = fopen(resname, 'r');
    res = fscanf(fid, '%d', Inf);
    fclose(fid);
    
    spkname = fullfile(spkfiles(i).folder, spkfiles(i).name);
    fid = fopen(spkname, 'r');
    spk = fread(fid, Inf, 'int16');
    spk = reshape(spk, nchans, nsamps, nspks);
    fclose(fid);
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % go over clusters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii = 1 : nclu
        
        % skip noise / artifact clusters
        if cluGrp(ii) == 0 || cluGrp(ii) == 1
            continue
        end
        
        % print progress
        if ii ~= 1
            fprintf(repmat('\b', 1, length(txt)))
        end
        txt = sprintf('Working on cluster %d / %d\n', ii, nclu);
        fprintf('%s', txt)
        
        % enforce dead time on spike detection
        cluIdx = find(clu == cluGrp(ii));
        nspksGrp = length(cluIdx);
        resGrp = res(cluIdx);  
        [resGrp, rmIdx, ~] = unique(resGrp);
        rmIdx = [setdiff(1 : nspksGrp, rmIdx)]';
        while min(diff(resGrp)) < dt
            ir = find(diff(resGrp) < dt);
            resGrp(ir) = [];
            rmIdx = [rmIdx; ir];
        end
        res(cluIdx(rmIdx)) = [];
        clu(cluIdx(rmIdx)) = [];
        spk(:, :, cluIdx(rmIdx)) = [];
        
        % recenter waveform peak
        cluIdx = find(clu == cluGrp(ii));
        nspksGrp = length(cluIdx);
        spkGrp = spk(:, :, cluIdx);
        spkMean = mean(spkGrp, 3);
        [~, ampMaxCh] = max(range(spkMean, 2));
        % here I should decide if inverted or regular spike and also think
        % of a way to only work on clusters where needed (e.g. big std) 
        
        for iii = 1 : nspksGrp
            [~, ia] = min(spkGrp(ampMaxCh, :, iii));
            iShift(iii) = psamp - ia;
            spkGrp(:, :, iii) = circshift(spkGrp(:, :, iii), iShift, 2);
        end
        
        % graphics
        figure
        subplot(2, 1, 1)
        stdshade(squeeze(spkGrp(ampMaxCh, :, :))', 0.3, 'k')
        subplot(2, 1, 2)
        stdshade(squeeze(spk(ampMaxCh, :, cluIdx))', 0.3, 'k')
        % also plot ccg 
        
        % update
        spk(:, :, cluIdx) = spkGrp;
        res(cluIdx) = res(cluIdx) + iShift;

                              
        % recalculate PCA
        
        
        % save files
        
        
    end
end

fprintf('that took %.2f minutes\n', toc / 60)





