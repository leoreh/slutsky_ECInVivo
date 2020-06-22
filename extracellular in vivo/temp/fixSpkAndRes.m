function fixSpkAndRes(varargin)

% processes neurosuite files after manual curation. answers two issues: (1)
% re-aligns individual waveforms according to their peak / trough
% (depending on the greater value in the mean waveform). (2) enforces a
% dead time on spike times (res) that typically occurs when kilosort tries
% to find spikes that occur simultaneously. re-alignment must occur first
% because it influences spike times. usually re-alignment also
% significantly improves the refractory period as seen in the ccg even
% without enforcing a dead time (increases confidance this is the right
% thing to do). alas recalculates pca. if backup is flagged then will copy
% the original files to a separate directory.
% -------------------------------------------------------------------------
% since the waveforms were detrended in ks2ns, we can assume the beginning
% and end of a wave are equal to zero and thus circhshift is used to
% recenter the peak / trough.
% -------------------------------------------------------------------------
% one persistant problem is that if a waveform includes two spikes, the
% larger one will always be centered even if the cluster originally
% referred to the smaller one.
% -------------------------------------------------------------------------
% 
% INPUT:
%   basepath    string. path to recording folder {pwd}. if multiple dat
%               files exist only the first will be processed
%   spkgrp      array where each cell is the electrodes for a spike group.
%               if empty will be loaded from session info file (cell
%               explorer format)
%   grp         numeric. groups (tetrodes) to work on. refers to filename
%               and not order. e.g. grp = 3 will load clu.3
%   fs          numeric. sampling frequency [hz]{20000}
%   dt          numeric. dead time to enforce between spikes [samples]{6}.
%   nsamps      numeric. length of waveform in spk file {32}. if empty will
%               be calculated according to fs as 1.6 ms 
%   psamp       numeric. peak / trough sample {16}. if empty will be
%               set to half nsamps.
%   stdFactor   numeric. re-alignment will occur only if std @ peak is
%               greater than std @ begining * stdFactor {1.1}. if empty
%               then all clusters will undergo re-alignment. 
%   graphics    logical. plot the before / after of each cluster {false}.
%   bkup        logical. create backup of files in basepath/ns/bkup {true}.
%   saveFiles   logical. save ns files after processing {true}. mainly for
%               debugging this code
%
% DEPENDENCIES
%   plotCCG
%   plotWaveform
%
% 22 jun 20 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'spkgrp', {}, @iscell);
addOptional(p, 'grp', [], @isnumeric);
addOptional(p, 'fs', 20000, @isnumeric);
addOptional(p, 'dt', 6, @isnumeric);
addOptional(p, 'nsamps', [], @isnumeric);
addOptional(p, 'psamp', [], @isnumeric);
addOptional(p, 'stdFactor', 1.1, @isnumeric);
addOptional(p, 'graphics', false, @islogical);
addOptional(p, 'bkup', true, @islogical);
addOptional(p, 'saveFiles', true, @islogical);

parse(p, varargin{:})
basepath    = p.Results.basepath;
spkgrp      = p.Results.spkgrp;
grp         = p.Results.grp;
fs          = p.Results.fs;
dt          = p.Results.dt;
nsamps      = p.Results.nsamps;
psamp       = p.Results.psamp;
stdFactor   = p.Results.stdFactor;
graphics    = p.Results.graphics;
bkup        = p.Results.bkup;
saveFiles   = p.Results.saveFiles;

% try to load params from session info (cell explorer)
[~, basename] = fileparts(basepath);
infoname = fullfile(basepath, [basename '.session.mat']);
if exist(infoname, 'file')
    load(infoname)
end
if isempty(spkgrp)
    spkgrp = session.extracellular.spikeGroups.channels;
end
ngrp = length(spkgrp);
if isempty(fs)
    fs = session.extracellular.sr;
end
if isempty(nsamps)
    nsamps = ceil(1.6 * 10^-3 * fs);
end
if isempty(psamp)
    psamp = nsamps / 2;
end
if isempty(grp)
    grp = 1 : ngrp;
end

% for CCG
dur = 0.01;
binsize = 0.0001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange files and backup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clufiles = dir([basepath filesep '*.clu*']);
clunames = {clufiles.name};
resfiles = dir([basepath filesep '*.res*']);
resnames = {resfiles.name};
spkfiles = dir([basepath filesep '*.spk*']);
spknames = {spkfiles.name};
fetfiles = dir([basepath filesep '*.fet*']);
fetnames = {fetfiles.name};
if length(resfiles) ~= length(clufiles) ||...
        length(clufiles) ~= length(spkfiles) ||...
        length(spkfiles) ~= length(fetfiles)
    error('different number of res / clu / spk / fet files')
end

% create backup
if bkup
    bkpath = fullfile(basepath, 'ns', 'bkup');
    fprintf('Saving backup in %s\n', bkpath)
    mkdir(bkpath)
    for i = grp
        copyfile(fullfile(clufiles(i).folder, clunames{i}), bkpath)
        copyfile(fullfile(resfiles(i).folder, resnames{i}), bkpath)
        copyfile(fullfile(spkfiles(i).folder, spknames{i}), bkpath)
        copyfile(fullfile(fetfiles(i).folder, fetnames{i}), bkpath)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% go over groups
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = grp
    nchans = length(spkgrp{i});
    
    % load data
    fprintf('Loading group %d / %d ', i, ngrp)
      
    cluname = fullfile(clufiles(i).folder, clufiles(i).name);
    fid = fopen(cluname, 'r');
    nclu = fscanf(fid, '%d', 1);
    clu = fscanf(fid, '%f', Inf);
    fclose(fid);
    cluGrp = unique(clu);
    if length(cluGrp) ~= nclu
        warning('cluFile has the wrong number of clusters')
        nclu = length(cluGrp);
    end
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
    txt = '';
    for ii = 1 : nclu
               
        % skip noise / artifact clusters
        if cluGrp(ii) == 0 || cluGrp(ii) == 1
            continue
        end

        % print progress
        fprintf(repmat('\b', 1, length(txt)))
        txt = sprintf('Working on cluster %d / %d \t', ii, nclu);
        fprintf('%s', txt)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % go over spikes and recenter 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cluIdx = find(clu == cluGrp(ii));  
        nspksGrp = length(cluIdx);
        spkGrp = spk(:, :, cluIdx);
        spkMean = mean(spkGrp, 3);
        [~, ampMaxCh] = max(range(spkMean, 2));
        
        % re-align waveforms only if std at peak greater than at edges
        spkStd = std(spkGrp(ampMaxCh, :, :), [], 3);
        if ~isempty(stdFactor)
            if spkStd(psamp) < spkStd(1) * stdFactor &&...
                    spkStd(psamp) < spkStd(end) * stdFactor
                continue
            end
        end
        
        % find wave polarity to decide if align to min or max
        if abs(min(spkMean(ampMaxCh, :))) > max(spkMean(ampMaxCh, :))
            spkPolarity = 'neg';
        else
            spkPolarity = 'pos';
        end
        
        % go over each spike and align to min / max
        if graphics
            spkOrig = spkGrp;
            resOrig = res(cluIdx);
        end
        iShift = zeros(nspksGrp, 1);
        for iii = 1 : nspksGrp
            switch spkPolarity
                case 'neg'
                    [~, ia] = min(spkGrp(ampMaxCh, :, iii));
                case 'pos'
                    [~, ia] = max(spkGrp(ampMaxCh, :, iii));
            end
            iShift(iii) = psamp - ia;
            spkGrp(:, :, iii) = circshift(spkGrp(:, :, iii), iShift(iii), 2);
        end       
        spk(:, :, cluIdx) = spkGrp;
        res(cluIdx) = res(cluIdx) + iShift;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % enforce dead time on spike detection
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        nspks = length(clu);

        % graphics
        if graphics
            figure
            txt2 = sprintf('Group %d Cluster %d', i, cluGrp(ii));
            suptitle(txt2)
            subplot(2, 2, 1)
            plotWaveform('avgwv', mean(spkOrig, 3),...
                'stdwv', std(spkOrig, [], 3), 'c', 'k',...
                'orient', 'vert', 'sbar', false)
            xticks([])
            yticks([])
            title('Before')
            subplot(2, 2, 2)
            plotWaveform('avgwv', mean(spkGrp(:, :, :), 3),...
                'stdwv', std(spkGrp(:, :, :), [], 3), 'c', 'k',...
                'orient', 'vert', 'sbar', false)
            xticks([])
            yticks([])
            title('After')
            subplot(2, 2, 3)
            [ccg, t] = CCG({resOrig / fs}, [], 'duration', dur, 'binSize', binsize);
            plotCCG('ccg', ccg(:, 1, 1), 't', t, 'basepath', basepath,...
                'saveFig', false, 'c', {'k'});
            subplot(2, 2, 4)
            cluIdx = clu == cluGrp(ii);
            [ccg, t] = CCG({res(cluIdx) / fs}, [], 'duration', dur, 'binSize', binsize);
            plotCCG('ccg', ccg(:, 1, 1), 't', t, 'basepath', basepath,...
                'saveFig', false, 'c', {'k'});
            drawnow;
        end                                        
    end
    
    % recalculate PCA
    fprintf('Calculating PCAs... \t')
    nFeatures = nchans * 3 + nchans + 1;
    fetMat = zeros(nspks, nFeatures);
    enrgIdx = nchans * 3;
    for ii = 1 : nchans
        [~, pcFeat] = pca(permute(spk(ii, :, :), [3, 2, 1]));
        chEnrgy = sum(abs(permute(spk(ii, :, :), [3, 2, 1])), 2);
        fetMat(:, ii * 3 - 2 : ii * 3) = (pcFeat(:, 1 : 3));
        fetMat(:, enrgIdx + ii) = (chEnrgy);
    end
    fetMat(:, end) = res;
    fet = int32(fetMat');
    
    % save files
    if saveFiles
        fprintf('Saving files... \t')
        fetname = fullfile(fetfiles(i).folder, fetnames{i});
        fid = fopen(fetname, 'w');
        formatstring = '%d';
        for ii = 2 : nFeatures
            formatstring = [formatstring, '\t%d'];
        end
        formatstring = [formatstring, '\n'];
        fprintf(fid, '%d\n', nFeatures);
        fprintf(fid, formatstring, fet);
        rc = fclose(fid);
        
        fid = fopen(spkname, 'w');
        fwrite(fid, spk(:), 'int16');
        rc = fclose(fid);
        
        fid = fopen(cluname, 'w');
        fprintf(fid, '%d\n', nclu);
        fprintf(fid, '%d\n', clu);
        rc = fclose(fid);
        
        fid = fopen(resname, 'w');
        fprintf(fid, '%d\n', res);
        rc = fclose(fid);
        fprintf('done\n')
    end
end

fprintf('that took %.2f minutes\n', toc / 60)

end

% EOF



