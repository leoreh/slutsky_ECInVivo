function numSpikes = getNumSpikes(basepath, spikes, saveVar)

% loads clu files and counts the number of spikes and clusters after
% automatic clustering and manual curation
%
% INPUT:
%   basepath        path to recording folder {pwd}.
%   saveVar         save variables
%
% OUTPUT:
%   numSpikes   struct
%
% 07 dec 18 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nargs = nargin;
if nargs < 1 || isempty(basepath)
    basepath = pwd;
end
if nargs < 2 || isempty(spikes)
    warning('spikes will be loaded from %s', basepath)
    spikes = getSpikes('basepath', basepath);
end
if nargs < 3 || isempty(saveVar)
    saveVar = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get clu files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(basepath)
clufiles = dir([basepath '\' '*.clu*']);
filenames = {clufiles.name};
nfiles = length(filenames);
if isempty(filenames)
    error('no .clu files in %s.', basepath)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and count spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : nfiles
    filename = fullfile(basepath, filenames{i});
    
    fid = fopen(filenames{i}, 'r');
    if(fid == -1);
        error('cannot open file');
    end
    nclu = fscanf(fid, '%d', 1);
    clu = fscanf(fid, '%f')';
    fclose(fid);
    
    % number of spikes per spike group
    numSpikes.spikesPostDetect(i) = length(clu);
    numSpikes.spikesPostCuration(i) = length(clu(clu > 1));
    numSpikes.ClustersPostCuration(i) = nclu - 2;
    
    if isfield(spikes, 'su')
        numSpikes.su(i) = sum(spikes.su' & spikes.shankID == i);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and count spikes from auto dir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist([basepath, '\auto'], 'dir')
    cd([basepath, '\auto'])   
    for i = 1 : nfiles
        filename = fullfile(basepath, filenames{i});
        
        fid = fopen(filenames{i}, 'r');
        if(fid == -1);
            error('cannot open file');
        end
        nclu = fscanf(fid, '%d', 1);
        clu = fscanf(fid, '%f')';
        fclose(fid);
        
        numSpikes.spikesPostAuto(i) = length(clu(clu > 1));
        numSpikes.ClustersPostAuto(i) = nclu - 2;

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saveVar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveVar
    cd(basepath)
    [~, basename] = fileparts(basepath);
    save([basename, '.numSpikes.mat'], 'clu');
end

end

% EOF