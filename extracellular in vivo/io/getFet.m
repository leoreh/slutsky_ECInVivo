function fet = getFet(basepath, saveVar, forceL)

% load features from all .fet files in folder.
%
% INPUT:
%   basepath        path to recording folder {pwd}.
%   saveVar         save output in basepath
%   forceL          reload fet even if .mat exists
%
% OUTPUT:
%   fet         array of k cells where k is the number of spike groups
%               (e.g. tetrodes). each cell is a matrix n x m + #elec + 1,
%               where n is the number of spikes and m is the number of
%               features. for each electrode there is one more column. Last
%               column is the clu ID array including noise (clu1) and
%               artifact (clu0) spikes.
%               
% 03 dec 18 LH
% 
% to do list:
%   divide last fet by fs
%   correct case were mat exists and forceL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nargs = nargin;
if nargs < 1 || isempty(basepath)
    basepath = pwd;
end
if nargs < 2 || isempty(saveVar)
    saveVar = 1;
end
if nargs < 3 || isempty(forceL)
    forceL = 0;
end

[~, filename, ~] = fileparts(basepath);
filename = [fullfile(basepath, filename) '.fet.mat'];
if exist(filename, 'file')  && ~forceL
    load(filename);
    warning('fet.mat already exists. loading %s', filename)
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get filenames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(basepath)
% .fet files
fetfiles = dir([basepath filesep '*.fet*']);
fetfilenames = {fetfiles.name};
nfetfiles = length(fetfilenames);
if isempty(fetfilenames)
    error('no .fet files in %s.', basepath)
end

% .clu files
clufiles = dir([basepath filesep '*.clu*']);
clufilenames = {clufiles.name};
nclufiles = length(clufilenames);
if isempty(clufilenames)
    error('no .fet files in %s.', basepath)
end

if nfetfiles ~= nclufiles
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
    fet{i}(:, end + 1) = fscanf(fid, '%f')';
    fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save fet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveVar
    save(filename, 'fet')
end

end

% EOF