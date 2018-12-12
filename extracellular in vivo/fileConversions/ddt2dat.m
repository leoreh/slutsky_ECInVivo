function ddt2dat(varargin)

% converts all ddt files (plexon stream) in basepath to dat files (neurosuite).
% performs basic preprocessing.
% moves ddt files to ddt folder.
% 
% INPUT:
%   basepath    path to recording folder {pwd}.
%   mapch       new order of channels {[]}.
%   rmvch       channels to remove (according to original order) {[]}
%   filenames   array of files to convert (including suffix .ddt)
%
% 22 nov 18 LH. updates:
% 30 nov 18 - handle multiple files
% 10 dec 18 - choose specific file
% 
% to do list:
% convert directly from TTank 
% option to concatenate files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'basepath', pwd, @ischar);
addOptional(p, 'mapch', []);
addOptional(p, 'rmvch', []);
addOptional(p, 'filenames', []);

parse(p, varargin{:});
mapch = p.Results.mapch;
rmvch = p.Results.rmvch;
basepath = p.Results.basepath;
filenames = p.Results.filenames;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get .ddt filenames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
cd(basepath)
if isempty(filenames)
    ddtfiles = dir([basepath '\' '*.ddt']);
    filenames = {ddtfiles.name};
    nfiles = length(filenames);
    if isempty(filenames)
        error('no .ddt files in %s.', basepath)
    end
else
    nfiles = size(filenames, 1);
end

fprintf(1, '\nFound %d .ddt files\n', nfiles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open files and extract info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : nfiles
    
    fprintf(1, 'Working on file %s\n', filenames{i});    
    [~, basename] = fileparts(filenames{i});
    newname = [basename '.dat'];
      
    fid = fopen(filenames{i}, 'r');
    fout = fopen(newname, 'w');  
    if(fid == -1);
        error('cannot open file');
    end
    if(fout == -1);
        error('cannot open file');
    end
    
    ddtversion = fread(fid, 1, 'int32');
    dataoffset = fread(fid, 1, 'int32');
    freq = fread(fid, 1, 'double');
    nch = fread(fid, 1, 'int32');
    fseek(fid, 0, 1);
    fsize = ftell(fid);
    frewind(fid);
    fseek(fid, dataoffset, 0);
    nsamples = (fsize - dataoffset)/(nch * 2);
    fclose(fid);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rearrange mapch according to rmvch
    % allows removing channels before remapping
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nchOrig = length(mapch);
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % partition into blocks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nbytes = 2;         % [bytes/sample]
    blocksize = 1e6;    % [samples/channel]
    
    % check file integrity
    if ~isequal(nsamples, round(nsamples))
        error('incorrect nchans (%d) for file %s', nch, filenames{i});
    end
    nblocks = ceil(nsamples / blocksize);
    blocks = [1 : blocksize : blocksize * nblocks; blocksize : blocksize : blocksize * nblocks]';
    blocks(nblocks, 2) = nsamples;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % go over blocks, remap, remove and write
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for j = 1 : nblocks
        boff = (blocks(j, 1) - 1) * nbytes * nch;
        bsize = (diff(blocks(j, :)) + 1);
        m = memmapfile(filenames{i}, 'Format', 'int16', 'Offset', boff, 'Repeat', bsize * nchOrig, 'writable', true);
        data = reshape(m.data, [nchOrig bsize]);
        
        if ~isempty(rmvch)                      % remove channels
            data(rmvch, :) = [];
        end
        
        if ~isempty(mapch)                      % remap channels
            data = data(mapch, :);
        end
        
        fwrite(fout, data(:), 'int16');         % write data
        
        clear data m
    end
    
    % close file
    rc = fclose(fout);
    if rc == -1
        error('cannot save new file');
    end
    
    info = dir(newname);
    fprintf(1, 'Created %s. File size = %.2f MB\n\n', newname, info.bytes / 1e6);
    
end

if ~exist('ddt', 'dir')
    sprintf('Creating ddt folder in %s', basepath)
    mkdir('ddt')
end
movefile('*.ddt', 'ddt');

fprintf(1, '\nElapsed time = %.2f seconds\n', toc)

end

% EOF