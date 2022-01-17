function datInfo = preprocDat(varargin)

% pre-process dat files. if no new desitation specified than creates backup
% and works on the new file. removes dc, removes pli, remaps and/or removes
% channels, handles acceleration, removes specified samples, creates info
% file. can concatenate multiple dat files.
%
% INPUT:
%   orig_files  char / cell of chars. full path and name of original dat
%               file/s. files will be concatenated based on their order in
%               filename
%   orig_paths  char / cell chars. full path where the original dat files
%               exist. will take all dat files in these paths and
%               concatenate them in the same order as they are naturally
%               sorted. if orig_files specified this will be ignored
%   newfile     char. full path and name of new file. if not
%               specified will be created in the folder of the first
%               datfile and will be named according to the datetime
%   chunksize   size of data to load at once [samples]{5e6}.
%               if empty will load entire file (be careful!).
%               for 35 channels in int16, 5e6 samples = 350 MB.
%   clip        cell of mat n x 2 indicating samples to diregard from chunks.
%               for example: clip{2} = [0 50; 700 Inf] will remove the first
%               50 samples and all samples between 700 and the end of the
%               2nd file.
%   nchans      numeric. original number of channels in the dat file {35}.
%   mapch       vec. new order of channels {[]}. 1-based.
%   rmvch       vec. channels to remove (according to original order) {[]}.
%               1-based.
%   pli         numeric. channel from which to extract line crossings.
%               if 0 will not remove pli.
%   precision   char. sample precision {'int16'}
%   saveVar     logical. save datInfo {true} or not (false).
%
% OUTPUT
%   datInfo     struct with fields describing original and processed files
%
% CALLS:
%   class2bytes
%   n2chunks
%
% TO DO LIST:
%   check linux compatible
%   conversion to mV
%   handle xml
%   datInfo (done)
%   add cat dat (done)
%   switch memmap for fread (done)
%
% 09 apr 20 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'orig_files', '');
addOptional(p, 'orig_paths', '');
addOptional(p, 'newfile', '', @ischar);
addOptional(p, 'chunksize', 5e6, @isnumeric);
addOptional(p, 'precision', 'int16', @ischar);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'clip', []);
addOptional(p, 'mapch', [], @isnumeric);
addOptional(p, 'rmvch', [], @isnumeric);
addOptional(p, 'pli', 0, @isnumeric);
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
orig_files = p.Results.orig_files;
orig_paths = p.Results.orig_paths;
newfile = p.Results.newfile;
chunksize = p.Results.chunksize;
precision = p.Results.precision;
nchans = p.Results.nchans;
clip = p.Results.clip;
mapch = p.Results.mapch;
rmvch = p.Results.rmvch;
pli = p.Results.pli;
saveVar = p.Results.saveVar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% size of one data point in bytes
nbytes = class2bytes(precision);

% check original files. if empty than get them from orig_paths
orig_files = cellstr(orig_files);
orig_paths = cellstr(orig_paths);
if isempty(orig_files{1})
    cnt = 1;
    for ipath = 1 : length(orig_paths)
        dat_files = dir([orig_paths{ipath} filesep '**' filesep '*dat']);
        for ifile = 1 : length(dat_files)
            orig_files{cnt} = fullfile(dat_files(ifile).folder, dat_files(ifile).name);
            cnt = cnt + 1;
        end
    end
end

% get info of original files
nfiles = length(orig_files);
for ifile = 1 : nfiles
    orig_info{ifile} = dir(orig_files{ifile});
    if isempty(orig_info{ifile})
        error('%s does not exist', orig_files{ifile})
    end
end

% handle clip
if isempty(clip)
    clip = cell(nfiles, 1);
end
if ~iscell(clip)
    clip = {clip};
end
if length(clip) ~= nfiles
    error('length of clip and orig_files not equal')
end
    
% create newname if doesn't exist
if isempty(newfile)
    [newpath, newname, ~] = fileparts(orig_files{1});
    newfile = fullfile(newpath, [newname, '_new.dat']);
end
[newpath, newname] = fileparts(newfile);

% rearrange mapch according to rmvch
orig_mapch = mapch;
if ~isempty(mapch)
    mapch = mapch(~ismember(mapch, rmvch));
    nch = length(mapch);
    if sum(mapch > nch)
        mat = [unique(mapch); (1 : nch)]';
        for j = 1 : nch
            mapch(j) = mat(mat(:, 1) == mapch(j), 2);
        end
    end
end
   
% open new file
fid = fopen(newfile, 'w');
if(fid == -1)
    error('cannot open file');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% go over chunks and process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ifile = 1 : nfiles
   
    fprintf('\npre-processing %s\n', orig_files{ifile})
 
    % open original file
    fid_in = fopen(orig_files{ifile}, 'r');
    if(fid_in == -1)
        error('cannot open file');
    end
    
    % partition into chunks
    nsamps(ifile) = orig_info{ifile}.bytes / nbytes / nchans;
    chunks = n2chunks('n', nsamps(ifile), 'chunksize', chunksize, 'clip', clip{ifile});
    nchunks = size(chunks, 1);

    % go over chunks
    for ichunk = 1 : nchunks
        
        % print progress
        if ichunk ~= 1
            fprintf(repmat('\b', 1, length(txt)))
        end
        txt = sprintf('working on chunk %d / %d', ichunk, nchunks);
        fprintf(txt)
       
        % load chunk
        dSamps = chunks(ichunk, 2) - chunks(ichunk, 1) + 1;
        dOffset = chunks(ichunk, 1) * nchans * nbytes - nbytes * nchans;
        fseek(fid_in, dOffset, 'bof');
        d = fread(fid_in, [nchans dSamps], [precision '=>' precision]);

        % remove channels
        if ~isempty(rmvch)
            % save channel used for line crossing detection
            if pli
                dpli = d(pli, :);
            end
            d(rmvch, :) = [];
        end
       
        % remap channels
        if ~isempty(mapch)
            d = d(mapch, :);
        end
       
        % remove dc
        [d, dc(ichunk, :)] = rmDC(d, 'dim', 2);
       
        % remove pli
        if pli
            linet = lineDetect('x', d(pli, :), 'fs', fs, 'graphics', false);
            for j = 1 : size(d, 1)
                d(j, :) = lineRemove(d(j, :), linet, [], [], 0, 1);
            end
        end
             
        % write to new file
        fwrite(fid, d(:), precision);
    end
   
    % close original file
    rc = fclose(fid_in);
    if rc == -1
        error('cannot close original file');
    end
end

% close new file
rc = fclose(fid);
fclose('all');
if rc == -1
    error('cannot save new file');
end

% check new file
info = dir(newfile);
nsampsNew = info.bytes / nbytes / (nchans - length(rmvch));
if ~isequal(nsampsNew, cumsum(nsamps))
    warning('processing failed, dats are of different length')
end
fprintf('\nfinished processing %s. \nFile size = %.2f MB\n', newfile, info.bytes / 1e6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange datInfo and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infoname = fullfile(newpath, [newname, '.datInfo.mat']);
if exist(infoname, 'file')
    load(infoname)
end

datInfo.orig_files = orig_files;
datInfo.orig_info = orig_info;
datInfo.nsamps = nsamps;
datInfo.clip = clip;
datInfo.mapCh = orig_mapch;
datInfo.rmvCh = rmvch;
datInfo.pli = pli;
datInfo.dc = mean(dc, 1);

if saveVar
    save(infoname, 'datInfo', '-v7.3');
end

fprintf('that took %.2f minutes\n', toc / 60)

end

% EOF