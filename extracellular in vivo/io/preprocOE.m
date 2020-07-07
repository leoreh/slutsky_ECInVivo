function datInfo = preprocOE(varargin)

% pre-process open ephys.
%
% INPUT:
%   basepath    string. path to recording folder {pwd}
%   exp         numeric. experiments to organize. if empty will process
%               all. refers to exp name and not number of folders (i.e. if
%               only exp 1 and 3 exist in basepath, than exp should be [1
%               3] and not [1 2].
%   rec         cell with n arrays where n = max(exp). each cell
%               contains indices to the relavent recordings of that
%               experiment. also corresponds to name and not number of
%               folders.
%   concat      logical. concatenate dat files in experiments {true}
%   mapch       vec. new order of channels {[]}. 1-based. 
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
%   # enable concatenation of experiments (done)
%   # add option to select specific recordings in experiments (done)
%   # fix exp idx (done)
%   # fix datenum (year sometimes wrong) (done 02 may 20)
%   # fix datenum in experiment name
%   # instead of saving OE tstamps, convert them to samples (Din)
%
% 09 apr 20 LH  updates
% 28 apr 20 LH  find exp by name and not number of folders
%               allowed user to select specific recordings
% 20 may 20 LH  concatenate experiments
%               copy xml

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'exp', [], @isnumeric);
addOptional(p, 'rec', {}, @iscell);
addOptional(p, 'concat', true, @islogical);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'mapch', [], @isnumeric);
addOptional(p, 'rmvch', [], @isnumeric);

parse(p, varargin{:})
basepath = p.Results.basepath;
exp = p.Results.exp;
rec = p.Results.rec;
concat = p.Results.concat;
nchans = p.Results.nchans;
mapch = p.Results.mapch;
rmvch = p.Results.rmvch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get files and folders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nprocessing %s\n\n', basepath)
cd(basepath)

% create new basefolder
[newpath, baseTime] = fileparts(basepath);
basename = bz_BasenameFromBasepath(newpath);
dn = datenum(baseTime, 'yyyy-MM-dd');
newpath = fullfile(newpath, [basename '_' datestr(dn, 'yyMMdd')]);
mkdir(newpath)

% find paths to relavent experiments
files = dir(['**' filesep '*.*']);
expFolders = files(contains({files.name}, 'experiment'));
expNames = natsort({expFolders.name});
if isempty(exp)      
    for i = 1 : length(expNames)
    exp(i) = str2num(char(regexp(expFolders(i).name, [filesep 'd*'], 'Match')));
    end
end
if isempty(rec)
    rec = cell(1, max(exp));
end

% find relavent recordings in experiment. 
k = 1;
for i = 1 : length(exp)
    expIdx(i) = find(strcmp(expNames, ['experiment' num2str(exp(i))]));
    exPath{i} = fullfile(expFolders(expIdx(i)).folder, expNames(expIdx(i)));
    % get time from xml file
    if isequal(exp(i), 1)
        xmlname = ['settings.xml'];
    else
        xmlname = ['settings_' num2str(exp(i)) '.xml'];
    end
    t = convert(xmltree(xmlname));
    dn = datenum(t.INFO.DATE, 'dd mmm yyyy HH:MM:SS');
    expTime{i} = datestr(dn, 'hhmmss');
    
    % get recording folders in each experiment
    recFolders = files(contains({files.name}, 'recording') &...
        strcmp({files.folder}, exPath{i}));
    recNames = natsort({recFolders.name});
    if isempty(rec{exp(i)})
        for iii = 1 : length(recNames)
            rec{exp(i)}(iii) = str2num(char(regexp(recNames{iii}, [filesep 'd*'], 'Match')));
        end
    end
    for ii = 1 : length(rec{exp(i)})
        recIdx(i) = find(strcmp(recNames, ['recording' num2str(rec{exp(i)}(ii))]));
        recPath{k} = fullfile(recFolders(recIdx(i)).folder, recFolders(recIdx(i)).name);
        k = k + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create new folder named according to time in xml file. assumes there
    % is one xml file per experiment named 'setting*'. the xml provides the
    % datetime of creation (i.e. recording) whereas file.date gets the
    % datetime of last edit
    expName = expTime{1};
    exPathNew = fullfile(newpath, expName);
    mkdir(exPathNew)
    fprintf('created %s\n', exPathNew)
    
    % move / concatenate dat
    catDat('basepath', recPath, 'newpath', exPathNew,...
        'concat', concat, 'nchans', nchans, 'saveVar', true);
    
    % pre-process dat
    datInfo = preprocDat('basepath', exPathNew, 'fname', '', 'mapch', mapch,...
        'rmvch', rmvch, 'nchans', nchans, 'saveVar', true,...
        'chunksize', 5e6, 'precision', 'int16', 'bkup', false);
    
    % get digital input
    getDinOE('basepath', recPath, 'newpath', exPathNew,...
        'concat', true, 'nchans', nchans, 'precision', 'int16',...
        'saveVar', true);
    
    % get acceleration
    newch = length(mapch) - length(rmvch);
    chAcc = [newch : -1 : newch - 2];
    EMGfromACC('basepath', exPathNew, 'fname', '',...
        'nchans', newch, 'ch', chAcc, 'force', false, 'saveVar', true,...
        'graphics', false, 'fsOut', 1250);
    
    % copy xml
    basefiles = dir(fileparts(newpath));
    xmlfiles = basefiles(contains({basefiles.name}, 'xml'));
    if length(xmlfiles) > 1
        for i = 1 : length(xmlfiles)
            if strcmp(xmlfiles(i).name, [basename, '.xml'])
                xmlfile = fullfile(xmlfiles(i).folder, xmlfiles(i).name);
                xmlnew = [exPathNew, filesep, expName '.xml'];
                copyfile(xmlfile, exPathNew)
                movefile([exPathNew, filesep, xmlfiles(i).name], xmlnew);
            end
        end
    elseif isempty(xmlfiles)
        fprintf('\nNo xml file in %s. skipping...\n',...
            fileparts(newpath))
    else
        xmlfile = fullfile(xmlfiles.folder, xmlfiles.name);
        xmlnew = [exPathNew, filesep, expName '.xml'];
        copyfile(xmlfile, exPathNew)
        movefile([exPathNew, filesep, xmlfiles.name], xmlnew);
    end
 
end


% EOF