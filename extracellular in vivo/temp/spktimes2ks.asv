function spktimes2ks(rez)

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
% 15 jun 20 LH.     updates:
% 07 aug 20 LH      fix case where no spikes in grp

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 20000;
peakSamp = [];
grps = [1 : 3];
basepath = pwd;
nchans = 35;

% constants   
sniplength = ceil(1.6 * 10^-3 * fs);
win = [-(floor(sniplength / 2) - 1) floor(sniplength / 2)];   
precision = 'int16'; % for dat file. size of one data point in bytes
nbytes = class2bytes(precision); 

if isempty(peakSamp)
    peakSamp = round(sniplength / 2);
end
if isempty(grps)
    grps = 1 : length(spkgrp);
end
ngrps = length(grps);

% build regressor for detrending
s = 0 : sniplength - 1;
scaleS = s(end);
a = s./scaleS;
b = max(a, 0);
W = b(:);
W = [reshape(W, sniplength, []), ones(sniplength,1)];
[Q, R] = qr(W,0);

% memory map to dat file
% get file
cd(basepath)
[~, basename] = fileparts(basepath);
fraw = dir([basename '.dat']);
nsamps = fraw.bytes / nbytes / nchans;
m = memmapfile(fraw.name, 'Format', {precision, [nchans, nsamps] 'mapped'});
raw = m.Data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% go over groups and save file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = 1;
for i = 1 : ngrps

    grp = grps(i);
    stamps{i} = spktimes{i}(:, 1);
    nspks(i) = length(stamps{i});
    grpchans = spkgrp{i};
    
    % ---------------------------------------------------------------------
    % spk file (binary, nsamples around each spike)
    spkname = fullfile([basename '.spk.' num2str(grp)]);
    fprintf('\nCreating \t%s. ', spkname)
    spk = zeros(length(grpchans), sniplength, nspks(i));   
    for ii = 1 : nspks(i)
        % print progress
        if mod(ii, 10000) == 0
            if ii ~= 10000
                fprintf(repmat('\b', 1, length(txt)))
            end
            txt = ['Extracted ', num2str(ii), ' / ', num2str(nspks(i)), ' spks'];
            fprintf('%s', txt)     
        end       
        % fix special case where spike is at end / begining of recording
        if stamps{i}(ii) + win(1) < 1 || stamps{i}(ii) + win(2) > nsamps
            warning('\nskipping stamp %d because waveform incomplete', ii)
            spk(:, :, ii) = [];
            stamps{i}(ii) = [];
            nspks(i) = nspks(i) - 1;
            continue
        end        
        % get waveform and remove best fit
        v = double(raw.mapped(grpchans, stamps{i}(ii) + win(1) :...
            stamps{i}(ii) + win(2)));
        v = [v' - W * (R \ Q' * v')]';

        % realign according to minimum
        [~, ib] = max(range(v'));   % this can be replaced by spktimes 2nd column, will be faster
        [~, ia] = min(v, [], 2);       
        peak = ia(ib);
        ishift = peak - peakSamp;
        if ishift ~= 0 
            k = k + 1;
            stamps{i}(ii) = stamps{i}(ii) + ishift;
                v = double(raw.mapped(grpchans, stamps{i}(ii) + win(1) :...
                    stamps{i}(ii) + win(2)));
                v = [v' - W * (R \ Q' * v')]';
        end
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
    % res 
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
    % clu (temporary, for visualization w/ neuroscope)
    nclu = 1;
    clugrp = ones(1, nspks(i));
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
    % fet file    
    fetname = fullfile([basename '.fet.' num2str(grp)]);
    fprintf('\nCreating \t%s. Computing PCAs...', fetname)
    nFeatures = length(grpchans) * 3 + length(grpchans) + 1;
    fetMat = zeros(nspks(i), nFeatures);
    enrgIdx = length(grpchans) * 3;
    if ~isempty(spk)
        for ii = 1 : length(grpchans)
            [~, pcFeat] = pca(permute(spk(ii, :, :), [3, 2, 1]));
            chEnrgy = sum(abs(permute(spk(ii, :, :), [3, 2, 1])), 2);
            fetMat(:, ii * 3 - 2 : ii * 3) = (pcFeat(:, 1 : 3));
            fetMat(:, enrgIdx + ii) = (chEnrgy);
        end
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