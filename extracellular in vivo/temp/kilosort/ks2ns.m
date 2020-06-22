function ks2ns(rez)

% converts KiloSort output (.rez structure) to neurosuite files (fet, res,
% clu, and spk). based in part on Kilosort2Neurosuite (Peterson). extracts
% and detrends waveforms from the dat file (see getSPKfromDat for
% more information). IMPROTANT: XML FILE MUST BE UPDATED TO THE PARAMS USED
% HERE (E.G. SAMPLES PER WAVEFORM)
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
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constants   
sniplength = ceil(1.6 * 10^-3 * rez.ops.fs);
win = [-(floor(sniplength / 2) - 1) floor(sniplength / 2)];   
precision = 'int16'; % for dat file. size of one data point in bytes
nbytes = class2bytes(precision); 

% params from rez
basepath = fileparts(rez.ops.fbinary);
[~, basename] = fileparts(basepath);
nchans = rez.ops.Nchan;
nchansTot = rez.ops.NchanTOT;

% spike times
if size(rez.st3, 2) > 4
    rez.st3 = rez.st3(:, 1 : 4);
end
[~, isort] = sort(rez.st3(:, 1), 'ascend');
rez.st3 = rez.st3(isort, :);
spktimes = uint64(rez.st3(:, 1));

% spike templates
rez.W = gather(single(rez.Wphy));
rez.U = gather(single(rez.U));
nt0 = size(rez.W, 1);
Nfilt = size(rez.W, 2);
templates = zeros(nchans, nt0, Nfilt, 'single');
for iNN = 1:size(templates, 3)
   templates(:, :, iNN) = squeeze(rez.U(:, iNN, :)) *...
       squeeze(rez.W(:, iNN, :))';
end
templates = permute(templates, [2, 1, 3]); % now it's nTemplates x nSamples x nChannels

% clusters
spikeTemplates = uint32(rez.st3(:, 2));
if size(rez.st3, 2) > 4
    spikeClusters = uint32(1 + rez.st3(:, 5));
    clu = uint32(spikeClusters - 1);
else
    clu = uint32(spikeTemplates - 1);
end

% channel map
chanmap = load(fullfile(basepath, 'ChanMap.mat'));
kc = chanmap.kcoords(chanmap.connected);

% find to which grp each template belongs
grps = unique(kc);
ngrps = length(grps);
for i = 1 : size(templates, 3)
    [~, ampMaxCh(i)] = max(range(templates(:, :, i)));
end

% build regressor for detrending
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
raw = m.Data;
clear rez isort spikeTemplates templates  % for memory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% go over groups and save file
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
            warning('\nskipping stamp %d because waveform incomplete', ii)
            spk(:, :, ii) = nan(length(chans), sniplength);
            continue
        end        
        % get waveform and remove best fit
        v = double(raw.mapped(chans, stamps{i}(ii) + win(1) : stamps{i}(ii) + win(2)));
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
     
end

% EOF