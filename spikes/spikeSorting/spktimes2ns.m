function spktimes2ns(varargin)

% creates ns files (spk, res, fet, and clu) from spktimes. can clip
% spktimes according to requested start time and/or duration. 
%
% INPUT:
%   basepath    string. path to recording folder {pwd}.
%   spkgrp      array where each cell is the electrodes for a spike group.
%               if empty will be loaded from session info file (cell
%               explorer format)
%   grps        numeric. groups (tetrodes) to work on
%   fs          numeric. sampling frequency [hz]{20000}
%   nchans      numeric. number of channels in dat file.
%   dur         numeric. duration of trim period [min]. 
%   t           string. start time of trim period. if empty than will take
%               dur minutes from end of recording (if dur is negative) or
%               dur minutes from start of recording (if dur is positive).
%               can be in the format 'HHmmss' or 'HHmm'.
%   mkClu       logical. create also clu file for inspection w/ ns {false}
%   spkFile     string. clip spike waveforms from 'dat' or {'temp_wh'}
%   saveVar     logical. save spktimes if updated during snipping. 
% 
% DEPENDENCIES
%   class2bytes
%   snipFromBinary
%
% TO DO LIST
%   # replace memmap w/ fread 
%   # clean temp_wh after creation of ns files
% 
% 09 nov 20 LH  updates
% 09 mar 20 LH  dur can be pos / neg
% 12 dec 21 LH  separated snip

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'spkgrp', {}, @iscell);
addOptional(p, 'grps', [], @isnumeric);
addOptional(p, 'fs', 20000, @isnumeric);
addOptional(p, 'nchans', [], @isnumeric);
addOptional(p, 'dur', [], @isnumeric);
addOptional(p, 't', []);
addOptional(p, 'mkClu', false, @islogical);
addOptional(p, 'spkFile', 'temp_wh', @ischar);
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
basepath    = p.Results.basepath;
spkgrp      = p.Results.spkgrp;
grps        = p.Results.grps;
fs          = p.Results.fs;
nchans      = p.Results.nchans;
dur         = p.Results.dur;
t           = p.Results.t;
mkClu   	= p.Results.mkClu;
spkFile   	= p.Results.spkFile;
saveVar   	= p.Results.saveVar;

if isempty(grps)
    grps = 1 : length(spkgrp);
end
ngrps = length(grps);

% handle binary file
[~, basename] = fileparts(basepath);
switch spkFile
    case 'dat'
        fname = fullfile(basepath, [basename, '.' spkFile]);
    case 'temp_wh'
        fname = fullfile(basepath, ['temp_wh.dat']);
        nchans = length([spkgrp{:}]);
end
info = dir(fname); 
nbytes = class2bytes('int16');
nsamps = info.bytes / nbytes / nchans;

% memory map to temp_wh binary
m = memmapfile(fname, 'Format', {'int16', [nchans, nsamps] 'mapped'});
raw = m.Data;

% snip params
sniplength = ceil(1.6 * 10^-3 * fs);
win = [-(floor(sniplength / 2) - 1) floor(sniplength / 2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trim spktimes according to dur and t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find boundry samples
dur = dur * 60 * fs;
recStart = '';
if ischar(t)
    recStart = split(basename, '_');
    recStart = recStart{end};
    if numel(recStart) == 6
        tformat = 'HHmmss';
    elseif numel(recStart) == 4
        tformat = 'HHmm';
    end
    recStart = datetime(recStart, 'InputFormat', tformat);
    t = datetime(t, 'InputFormat', tformat);
    if t <= recStart
        t = t + hours(24);
    end
    s = seconds(t - recStart) * fs;
    if ~isempty(dur)
        trimEdges = sort([s s + dur]);
    else
        trimEdges = [s nsamps];
    end
else
    if dur < 0 
        trimEdges = [nsamps + dur nsamps];
    elseif dur > 0
        trimEdges = [1 dur];
    else
        trimEdges = [1 nsamps];
    end
end
if trimEdges(1) < 1
    warning('Requested time beyond recording duration')
    trimEdges(1) = 1;
elseif trimEdges(2) > nsamps
    warning('Requested time beyond recording duration')
    trimEdges(2) = nsamps;
end

% trim spktimes
load([basename '.spktimes.mat'])
for igrp = 1 : ngrps
    spktimes{igrp} = spktimes{igrp}(spktimes{igrp} > trimEdges(1)...
        & spktimes{igrp} < trimEdges(2));
end

% save trim info 
infoname = fullfile(basepath, [basename, '.datInfo.mat']);
if exist(infoname, 'file')
    load(infoname)
    datInfo.spktrim.dur = dur / fs / 60;
    datInfo.spktrim.edges = trimEdges;
    datInfo.spktrim.recStart = recStart;
    datInfo.spktrim.t = t;
    save(infoname, 'datInfo', '-v7.3');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% go over groups and save file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for igrp = 1 : ngrps
    
    grp = grps(igrp);
    grpchans = spkgrp{grp};

    % snip spikes from whitened data
    [spks, spktimes{grp}] = snipFromBinary('stamps', spktimes{grp}, 'fname', '',...
        'win', win, 'nchans', nchans, 'ch', grpchans, 'align_peak', 'min',...
        'precision', 'int16', 'rmv_trend', 0, 'saveVar', false, 'raw', raw);
        
    % ---------------------------------------------------------------------
    % spk file (binary, nsamples around each spike)   
    nspks = length(spktimes{grp});
    spkname = fullfile([basename '.spk.' num2str(grp)]);
    fid = fopen(spkname, 'w');
    fwrite(fid, spks(:), 'int16');
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
    fprintf(fid, '%d\n', spktimes{grp});
    rc = fclose(fid);
    if rc == 0
        fprintf('\nCreated \t%s', resname)
    else
        fprintf('\nFailed to create %s!', resname)
    end

    % ---------------------------------------------------------------------
    % clu (temporary, for visualization w/ neuroscope)
    if mkClu
        mkdir(['kk' filesep 'preSorting']) 
        nclu = 1;
        clugrp = ones(1, nspks);
        cluname = fullfile(['kk' filesep 'preSorting'], [basename '.clu.' num2str(grp)]);
        fid = fopen(cluname, 'w');
        fprintf(fid, '%d\n', nclu);
        fprintf(fid, '%d\n', clugrp);
        rc = fclose(fid);
        if rc == 0
            fprintf('\nCreated \t%s', cluname)
        else
            fprintf('\nFailed to create %s!', cluname)
        end
    end
    
    % ---------------------------------------------------------------------
    % fet file    
    fetname = fullfile([basename '.fet.' num2str(grp)]);
    fprintf('\nCreating \t%s. Computing PCAs...', fetname)
    nFeatures = length(grpchans) * 3 + length(grpchans) + 1;
    fetMat = zeros(nspks, nFeatures);
    enrgIdx = length(grpchans) * 3;
    if ~isempty(spks)
        for ichan = 1 : length(grpchans)
            [~, pcFeat] = pca(permute(spks(ichan, :, :), [3, 2, 1]), 'NumComponents', 3);
            chEnrgy = sum(abs(permute(spks(ichan, :, :), [3, 2, 1])), 2);
            fetMat(:, ichan * 3 - 2 : ichan * 3) = pcFeat;
            fetMat(:, enrgIdx + ichan) = (chEnrgy);
        end
    end
    fetMat(:, end) = double(spktimes{grp});
    fet = int32(fetMat');
    fid = fopen(fetname, 'w');
    formatstring = '%d';
    for ifet = 2 : nFeatures
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

% save variable
if saveVar
    save(fullfile(basepath, [basename, '.spktimes.mat']), 'spktimes')
end

clear raw
clear m

end

% EOF