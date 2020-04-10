function preprocOE(varargin)

% pre-process open ephys. 
%
% INPUT:
%   basepath    string. path to recording folder {pwd}
%   exp         numeric. experiments to organize. if empty will process all
%   concat      logical. concatenate dat files in experiments {true}
%   mapch       vec. new order of channels {[]}
%   rmvch       vec. channels to remove (according to original order) {[]}
%   nchans      numeric. number of channels in dat file {35}.
%
% CALLS:
%   xmltree
%   catDat
%   preprocDat
%   getDinOE
%
% TO DO LIST:
%   enable concatenation of experiments (not only recordings)
%
% 09 apr 20 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'exp', [], @isnumeric);
addOptional(p, 'concat', true, @islogical);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'mapch', [], @isnumeric);
addOptional(p, 'rmvch', [], @isnumeric);

parse(p, varargin{:})
basepath = p.Results.basepath;
exp = p.Results.exp;
concat = p.Results.concat;
nchans = p.Results.nchans;
mapch = p.Results.mapch;
rmvch = p.Results.rmvch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get files and folders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(basepath)

% create new basefolder
[newpath, baseTime] = fileparts(basepath);
basename = bz_BasenameFromBasepath(newpath);
dn = datenum(baseTime, 'yyyy-MM-dd_hh-mm-ss');
newpath = fullfile(newpath, [basename '_' datestr(dn, 'yyMMdd')]);
mkdir(newpath)

% get files
files = dir(['**' filesep '*.*']);
datFiles = files(contains({files.name}, '.dat'));
expFolders = files(contains({files.name}, 'experiment'));
if isempty(exp)
    exp = 1 : length(expFolders);
end
fprintf('\nfound %d experiments\n', length(expFolders))

% rename experiments based on time in xml file. assumes there is one xml
% file per experiment named 'setting*'. the xml provides the datetime of
% creation (i.e. recording) whereas file.date gets the datetime of last
% edit
for i = exp
    
    % print to 
    fprintf('working on experiment %d\n', i)
    
    % exp path and number
    expPath = fullfile(expFolders(i).folder, expFolders(i).name);
    expIdx = char(regexp(expFolders(i).name, [filesep 'd*'], 'Match'));
    
    % get time from xml file
    if strcmp(expIdx, '1')
        xmlname = ['settings.xml'];
    else
        xmlname = ['settings_' expIdx '.xml'];
    end
    t = convert(xmltree(xmlname));
    expTime = datestr(t.INFO.DATE, 'hhmmss');
    
    % get recording folders in each experiment 
    recFolders = files(contains({files.name}, 'recording') &...
        strcmp({files.folder}, expPath));
    recIdx = [];
    for j = 1 : length(recFolders)
            recIdx(j) = str2num(char((regexp(recFolders(j).name, [filesep 'd*'], 'Match'))));
    end
    
    fprintf('%d recordings found\n', length(recIdx))
    
    if concat
        expName = sprintf('%s_e%sr%d-%d', expTime, expIdx, recIdx(1), recIdx(end));
        expPathNew = fullfile(newpath, expName);
        mkdir(expPathNew)        
        fprintf('created %s\n', expPathNew)
        
        catDat('datpath', expPath, 'newpath', expPathNew,...
            'concat', true, 'nchans', 35, 'saveVar', false);       
        preprocDat('datpath', expPathNew, 'fname', '', 'mapch', [],...
            'rmvch', [1 : 20], 'nchans', 35, 'saveVar', false,...
            'chunksize', 5e6, 'precision', 'int16', 'bkup', true);
        getDinOE('datpath', expPath, 'newpath', expPathNew,...
            'concat', true, 'nchans', 35, 'nbytes', 2,...
            'saveVar', true);
    else
        for j = 1 : length(recFolders)
            % make new dir
            expName = sprintf(' %s_e%s_r%s', expTime, expIdx, recIdx(j));
            expPathNew = fullfile(newpath, expName);
            mkdir(expPathNew)
            % get dat
            recPath = fullfile(recFolders(j).folder, recFolders(j).name);
            preprocDat('basepath', recPath, 'newpath', expPathNew, 'concat', false)
            % Din
        end
    end
end

end

    
    

% EOF