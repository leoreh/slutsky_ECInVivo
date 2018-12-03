function fet = getFet(basepath)

% load features from all .fet files in folder.
%
% INPUT:
%   basepath    path to recording folder {pwd}.
%   saveMat     save output in basepath
%
% OUTPUT:
%   fet         cell array of k x n x m + 1, where k is the number of spike
%               groups (e.g. tetrodes), n is the number of spikes and m is
%               the number of features. The last column is the clu ID
%               
% 03 dec 18 LH
% 
% to do list:
%   divide last fet by fs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nargs = nargin;
if nargs < 1 || isempty(basepath)
    basepath = pws;
end
if nargs < 2 || isempty(saveMat)
    saveMat = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get filenames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(basepath)
% .fet files
fetfiles = dir([basepath '\' '*.fet*']);
fetfilenames = {fetfiles.name};
nfetfiles = length(fetfilenames);
if isempty(fetfilenames)
    error('no .fet files in %s.', basepath)
end

% .clu files
clufiles = dir([basepath '\' '*.clu*']);
clufilenames = {clufiles.name};
nclufiles = length(clufilenames);
if isempty(clufilenames)
    error('no .fet files in %s.', basepath)
end

if length(nfetfiles) ~= length(nclufiles)
    error('different number of .fet (%d) and clu (%d) files.',...
        length(nfetfiles), length(nclufiles))
end

fprintf(1, '\nFound %d .fet & .clu files\n', nclufiles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load fet and clu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : nclufiles
    
    fprintf(1, 'Working on file %s\n', fetfilenames{i});
    
    fid = fopen(fetfilenames{i}, 'r');
    if(fid == -1);
        error('cannot open file');
    end
    nfet = fscanf(fid, '%d', 1);
    fet{i} = fscanf(fid, '%f', [nfet, inf])';
    fclose(fid);
      
    fid = fopen(clufilenames{i}, 'r');
    if(fid == -1);
        error('cannot open file');
    end
    nclu = fscanf(fid, '%d', 1);
    fet{i}(:, end+1) = fscanf(fid, '%f')';
    fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save fet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveMat
    [~, filename, ~] = fileparts(basepath);
    save([fullfile(basepath, filename) '.fet.mat'], 'fet')
end

end

% EOF