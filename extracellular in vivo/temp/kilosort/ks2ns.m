function ks2ns(rez)

% converts KiloSort output (.rez structure) to neurosuite files (fet, res,
% clu, and spk). based in part on Kilosort2Neurosuite (Peterson). extracts
% and detrends waveforms from the dat file (see getSPKfromDat for
% moreinformation). calculated pca via parallel computing. loads all
% waveforms to memory at once. note that in order for this to work the xml
% file must contain spike groups with the correct number of samples. 
%
% DEPENDENCIES
%   class2bytes
%
% TO DO LIST:
%   # find trend based on larger segment without spike
%   # fix clipped spike at end / beginning of recording
%
% 15 june 20 LH

t1 = tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constants
win = [-19 20]; % win for snipping spk    
precision = 'int16'; % for dat file. size of one data point in bytes
nbytes = class2bytes(precision); 

% from rez
basepath = fileparts(rez.ops.fbinary);
[~, basename] = fileparts(basepath);
nchans = rez.ops.Nchan;
nchansTot = rez.ops.NchanTOT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load npy files and channel map
clu = readNPY(fullfile(basepath, 'spike_clusters.npy'));
spktimes = readNPY(fullfile(basepath, 'spike_times.npy'));
templates = readNPY(fullfile(basepath, 'templates.npy'));
templates = permute(templates,[2,3,1]);
chanmap = load(fullfile(basepath, 'ChanMap.mat'));
kc = chanmap.kcoords(chanmap.connected);
% spktimes        = uint64(rez.st3(:, 1));
% spktemplates    = uint32(rez.st3(:, 2)); % template id for each spike

% find to which grp each template belongs
grps = unique(kc);
ngrps = length(grps);
for i = 1 : size(templates, 3)
    [~, ampMaxCh(i)] = max(range(templates(:, :, i)));
end

% build regressor for detrending
sniplength = length(win(1) : win(2));
s = 0 : sniplength - 1;
scaleS = s(end);
a = s./scaleS;
b = max(a, 0);
W = b(:);
W = [reshape(W, sniplength, []), ones(sniplength,1)];
[Q, R] = qr(W,0);

% memory map to dat file
datname = rez.ops.fbinary;
info = dir(datname);
nsamps = info.bytes / nbytes / nchansTot;
m = memmapfile(datname, 'Format', {precision, [nchansTot, nsamps] 'mapped'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save files for each group
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : ngrps
    grp = grps(i);
    inShank = ismember(ampMaxCh, find(kc == i));
    shankClu = find(inShank == 1) - 1;
    idx = ismember(clu, shankClu);
    clugrp = clu(idx) + 2;
    nclu = length(shankClu);
    chans = find(kc == i);
    stamps{i} = spktimes(idx);
    nspks(i) = length(stamps{i});
    
    % ---------------------------------------------------------------------
    % res (ASCII file, one row per spike time [samples])
    resname = fullfile([basename '.res.' num2str(grp)]);
    fid = fopen(resname, 'w');
    fprintf(fid, '%d\n', stamps{i});
    rc = fclose(fid);
    if rc == 0
        fprintf('\nCreated \t%s', resname)
    else
        fprintf('\nFailed to create %s!', resname)
    end

    % ---------------------------------------------------------------------
    % clu file (ASCII file, one row per clu ID. first row is number of
    % clusters)
    cluname = fullfile([basename '.clu.' num2str(grp)]);
    fid = fopen(cluname, 'w');
    fprintf(fid, '%d\n', nclu);
    fprintf(fid, '%d\n', clugrp);
    rc = fclose(fid);
    if rc == 0
        fprintf('\nCreated \t%s', cluname)
    else
        fprintf('\nFailed to create %s!', cluname)
    end
    
    % ---------------------------------------------------------------------
    % spk file (binary, nsamples around each spike)
    spkname = fullfile([basename '.spk.' num2str(grp)]);
    fprintf('\nCreating \t%s. ', spkname)
    spk = zeros(length(chans), sniplength, nspks(i));   
    for ii = 1 : nspks(i)
        % print progress
        if mod(ii, 10000) == 0
            if ii ~= 10000
                fprintf(repmat('\b', 1, length(txt)))
            end
            txt = ['Extracted ', num2str(ii), ' / ', num2str(nspks(i)), ' spks'];
            fprintf('%s', txt)     
        end
        
        % lazy fix to special case where spike is at end of recording. will
        % not work because must update clu and res file.
        if stamps{i}(ii) + win(1) < 1 || stamps{i}(ii) + win(2) > nsamps
            warning('skipping stamp %d because waveform incomplete', i)
            spk(:, :, i) = nan(length(chans), sniplength);
            continue
        end
        
        % get waveform and remove best fit
        v = double(m.Data.mapped(chans,...
            stamps{i}(ii) + win(1) : stamps{i}(ii) + win(2)));
        v = [v' - W * (R \ Q' * v')]';
        spk(:, :, ii) = v;
    end
    % save to spk file
    fid = fopen(spkname, 'w');
    fwrite(fid, spk(:), 'int16');
    rc = fclose(fid);
    if rc == 0
        fprintf('. done')
    else
        fprintf('. Failed to create %s!', spkname)
    end
    
    % ---------------------------------------------------------------------
    % fet file    
    fetname = fullfile([basename '.fet.' num2str(grp)]);
    fprintf('\nCreating \t%s. Computing PCAs...', fetname)
    nFeatures = length(chans) * 3 + length(chans) + 1;
    fetMat = zeros(nspks(i), nFeatures);
    enrgIdx = length(chans) * 3;
    for ii = 1 : length(chans)
       [~, pcFeat] = pca(permute(spk(ii, :, :), [3, 2, 1]));
       chEnrgy = sum(abs(permute(spk(ii, :, :), [3, 2, 1])), 2);
       fetMat(:, ii * 3 - 2 : ii * 3) = (pcFeat(:, 1 : 3));
       fetMat(:, enrgIdx + ii) = (chEnrgy);
    end
    fetMat(:, end) = double(stamps{i});
    fet = int32(fetMat');
    fid = fopen(fetname, 'w');
    formatstring = '%d';
    for ii = 2 : nFeatures
        formatstring = [formatstring, '\t%d'];
    end
    formatstring = [formatstring, '\n'];  
    fprintf(fid, '%d\n', nFeatures);
    fprintf(fid, formatstring, fet);
    rc = fclose(fid);
    if rc == 0
        fprintf(' done\n')
    else
        fprintf(' Failed to create %s\n!', fetname)
    end
end
     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % res and clu
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % assign units (templates) to groups
% templates = gpuArray(zeros(nchans, size(rez.W, 1), size(rez.U, 2), 'single'));
% for i = 1 : size(rez.U, 2)
%     templates(:, :, i) = squeeze(rez.U(:, i, :)) * squeeze(rez.W(:, i, :))';
% end
% for i = 1 : size(templates, 3)
%     [~, ampMaxCh(i)] = max(range(templates(:, :, i)'));
% end
% 
% unitGrp = kcoords(ampMaxCh);
% grps = unique(unitGrp);
% ngrp = length(grps);
% 
% for i = 1 : ngrp
%     grp = grps(i);
%     
%     % res files
%     inShank          = ismember(ampMaxCh, find(kcoords == i));
%     shankClu         = find(inShank == 1) - 1;
%     idx              = ismember(clu, shankClu);
%     tclu             = clu(idx) + 2;
%     stamps{i}        = res(idx);
%     
%     
%     fprintf(['\nSaving .res and .clu file for group ', num2str(grp)])
%     resname = fullfile([basename '.res.' num2str(grp)]);
%     templateIdx = find(unitGrp == grp);
%     ia{i} = find(ismember(spktemplates, templateIdx));
%     stamps{i} = spktimes(ia{i});      % spk times per group
%     fid = fopen(resname, 'w');
%     fprintf(fid, '%d\n', stamps{i});
%     fclose(fid);
%     clear fid
%     
%     % clu files
%     cluname = fullfile([basename '.clu.' num2str(grp)]);
%     tclu = spktemplates(ia{i});
%     fid = fopen(cluname, 'w');
%     fprintf(fid, '%d\n', length(unique(tclu)));
%     fprintf(fid, '%d\n' , uint32(tclu));
%     fclose(fid);
%     clear fid
% end
% fprintf('\n'); toc(t1)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % spk
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% % snip waveforms and detrend
% for i = 1 : ngrp
%     grp = grps(i);
%     fprintf('\nExtracting waveforms for group %d', grp)    
%     chans = find(grp == kcoords);
%     nspks(i) = length(stamps{i});
%     wv{i} = zeros(length(chans), sniplength, nspks(i));
%     
%     for ii = 1 : nspks(i)
%         if mod(ii, 10000) == 0
%             if ii ~= 10000
%                 fprintf(repmat('\b', 1, length(txt)))
%             else
%                 fprintf('\n')
%             end
%             txt = ['finished ', num2str(ii), ' / ', num2str(nspks(i))];
%             fprintf('%s', txt)     
%         end
%         
%         % lazy fix to special case where spike is at end of recording. will
%         % not work because must update clu and res file.
%         if stamps{i}(ii) + win(1) < 1 || stamps{i}(ii) + win(2) > nsamps
%             error('skipping stamp %d because snip incomplete', i)
%             snips(:, :, i) = nan(length(chans), sniplength);
%             continue
%         end
%         
%         % get waveform and remove best fit
%         v = double(m.Data.mapped(chans,...
%             stamps{i}(ii) + win(1) : stamps{i}(ii) + win(2)));
%         v = [v' - W * (R \ Q' * v')]';
%         
%         wv{i}(:, :, ii) = v;
%     end
%     
%     % save to spk file
%     spkname = fullfile([basename '.spk.' num2str(grp)]);
%     fid = fopen(spkname, 'w');
%     fwrite(fid, wv{i}(:), 'int16');
%     fclose(fid);
%     clear fid
% end
% fprintf('\n'); toc(t1)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fet
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % start parpool
% % if ~parcompute && isempty(gcp('nocreate'))
% %     parpool
% % end
% 
% for i = 1 : ngrp
%     
%     spk = wv{i};
%     chans = chans;
%     
%     nFeatures = length(chans) * 3 +length(chans) + 1;
%     fetMat = zeros(nspks(1), nFeatures);
%     enrgIdx = length(chans) *3 ;
%     for i = 1 : length(chans)
%        [~, pcFeat] = pca(permute(spk(i, :, :), [3, 2, 1]));
%        chEnrgy = sum(abs(permute(spk(i, :, :), [3, 2, 1])), 2);
%        fetMat(:, i * 3 - 2 : i * 3) = (pcFeat(:, 1 : 3));
%        fetMat(:, enrgIdx + i) = (chEnrgy);
%     end
%     fetMat(:, end) = double(stamps{i});
%     fet                         = int32(fetMat');
%     fid=fopen(fetfname,'w');
%     formatstring = '%d';
%     for ii=2:nFeatures
%         formatstring = [formatstring,'\t%d'];
%     end
%     formatstring = [formatstring,'\n'];
%     
%     fprintf(fid, '%d\n', nFeatures);
%     fprintf(fid,formatstring,fet);
%     fclose(fid);
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     grp = grps(i);
%     wvforms = wv{i};   
%     nFeatures = length(chans) * 3 + length(chans) + 1;
%     fetMat = zeros(nspks(i), nFeatures);
%     enrgIdx = length(chans) * 3;
%     for ii = 1 : length(chans)
%        [~, pcFeat] = pca(permute(wvforms(ii, :, :), [3, 2, 1]));
%        chEnrgy = sum(abs(permute(wvforms(ii, :, :), [3, 2, 1])), 2);
%        fetMat(:, ii * 3 - 2 : ii * 3) = (pcFeat(:, 1 : 3));
%        fetMat(:,enrgIdx + ii) = (chEnrgy);
%     end
%     fetMat(:, end) = double(stamps{i});
%     fet = int32(fetMat');
%     
%     % save to fet
%     fetname = fullfile([basename '.fet.' num2str(grp)]);
%     fid = fopen(fetname, 'w');
%     formatstring = '%d';
%     for ii = 2 : nFeatures
%         formatstring = [formatstring, '\t%d'];
%     end
%     formatstring = [formatstring, '\n'];   
%     fprintf(fid, '%d\n', nFeatures);
%     fprintf(fid, formatstring, fet);
%     fclose(fid);
%     
%   
% %     grp = grps(i);
% %     fprintf(['\nComputing PCAs for group ', num2str(grp),])
% %     pcaTotal = zeros(3, length(chans), length(stamps{i}));
% %     chEnergy = zeros(length(chans), length(stamps{i}));
% %     wvforms = wv{i};
% %     
% %     % calculate pca
% %     if isempty(gcp('nocreate'))
% %         for ii = 1 : length(chans)
% %             [~, pcaCh] = pca(zscore(permute(wvforms(ii, :, :),...
% %                 [3, 2, 1]), [], 2), 'NumComponents', 3);
% %             pcaTotal(:, ii, :) = pcaCh';
% %             chEnergy(ii, :) = sum(abs(permute(wvforms(ii, :, :), [3, 2, 1])), 2);
% %         end
% %     else
% %         parfor ii = 1 : length(chans)
% %             [~, pcaCh] = pca(zscore(permute(wvforms(ii, :, :),...
% %                 [3, 2, 1]), [], 2), 'NumComponents', 3);
% %             pcaTotal(:, ii, :) = pcaCh';
% %             chEnergy(ii, :) = sum(abs(permute(wvforms(ii, :, :), [3, 2, 1])), 2);
% %         end
% %     end
% %     pcaTotal = reshape(pcaTotal, size(pcaTotal, 1) *...
% %         size(pcaTotal, 2), size(pcaTotal, 3));
% %     
% %     fet = int32([pcaTotal; chEnergy; stamps{i}']);
% %     nFeatures = size(fet, 1);
% %     
% %     % save to fet
% %     formatstring = '%d';
% %     for ii = 2 : nFeatures
% %         formatstring = [formatstring, '\t%d'];
% %     end
% %     formatstring = [formatstring, '\n'];
% %     fetname = fullfile([basename '.fet.' num2str(grp)]);
% %     fid = fopen(fetname, 'w');
% %     fprintf(fid, '%d\n', nFeatures);
% %     fprintf(fid, formatstring, fet);
% %     fclose(fid);
% end
% 
% fprintf('\n'); toc(t1)
% fprintf('\ndone!\n')

end

% EOF