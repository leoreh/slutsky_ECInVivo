function datInfo = preprocOE(varargin)

% pre-process open ephys.
%
% INPUT:
%   basepath    string. path to recording folder {pwd}
%   exp         numeric. experiments to organize. if empty will process all
%               refers to exp name and not number of folders (i.e. if only
%               exp 1 and 3 exist in basepath, than exp should be [1 3] and
%               not [1 2].
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
%   # enable concatenation of experiments (not only recordings)
%   # add option to select specific recordings in experiments
%   # fix exp idx (done)
%   # fix datenum (year sometimes wrong) (done 02 may 20)
%
% 09 apr 20 LH  updates
% 28 apr 20 LH  find exp by name and not number of folders
%               allowed user to select specific recordings

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

fprintf('\npre-processing %s\n\n', basepath)
cd(basepath)

% create new basefolder
[newpath, baseTime] = fileparts(basepath);
basename = bz_BasenameFromBasepath(newpath);
dn = datenum(baseTime, 'yyyy-MM-dd');

newpath = fullfile(newpath, [basename '_' datestr(dn, 'yyMMdd')]);
mkdir(newpath)

% find experiments
files = dir(['**' filesep '*.*']);
expFolders = files(contains({files.name}, 'experiment'));
if isempty(exp)
    exp = 1 : length(expFolders);
end
fprintf('\nfound %d experiments\n', length(expFolders))

% rename experiments based on time in xml file. assumes there is one xml
% file per experiment named 'setting*'. the xml provides the datetime of
% creation (i.e. recording) whereas file.date gets the datetime of last
% edit
for i = 1 : length(exp)
    
    fprintf('working on experiment %d\n', exp(i))
    
    % exp path
    expNames = {expFolders.name};
    expIdx = find(strcmp(expNames, ['experiment' num2str(exp(i))]));
    exPath = fullfile(expFolders(expIdx).folder, expFolders(expIdx).name);
    expIdx = char(regexp(expFolders(expIdx).name, [filesep 'd*'], 'Match'));
    
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
        strcmp({files.folder}, exPath));
    recIdx = [];
    for j = 1 : length(recFolders)
        recIdx(j) = str2num(char((regexp(recFolders(j).name, [filesep 'd*'], 'Match'))));
    end
    
    fprintf('%d recordings found\n', length(recIdx))
    
    expName = sprintf('%s_e%sr%d-%d', expTime, expIdx, recIdx(1), recIdx(end));
    exPathNew = fullfile(newpath, expName);
    mkdir(exPathNew)
    fprintf('created %s\n', exPathNew)
    
    % move / concatenate dat
    catDat('basepath', exPath, 'newpath', exPathNew,...
        'concat', concat, 'nchans', nchans, 'saveVar', true);
    
    % pre-process dat
    datInfo = preprocDat('basepath', exPathNew, 'fname', '', 'mapch', mapch,...
        'rmvch', rmvch, 'nchans', nchans, 'saveVar', true,...
        'chunksize', 5e6, 'precision', 'int16', 'bkup', false);
    
    % get digital input
    getDinOE('datpath', exPath, 'newpath', exPathNew,...
        'concat', true, 'nchans', nchans, 'nbytes', 2,...
        'saveVar', true);
    
    % get acceleration
    newch = length(mapch) - length(rmvch);
    chAcc = [newch : -1 : newch - 2];
    getAccFromDat('basepath', exPathNew, 'fname', '',...
        'nchans', newch, 'ch', chAcc, 'force', false, 'saveVar', true,...
        'graphics', false, 'fs', 1250);    
end

end

% EOF