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
%   mapch       vec. new order of channels {[]}. 1-based. 
%   rmvch       vec. channels to remove (according to original order) {[]}
%   nchans      numeric. number of channels in dat file {35}.
%   accCh       vec. acceleration channels. must be according to new channel 
%               order. if empty will skip {[]}.
%   fsIn        numeric. sampling frequency of recording
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
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'mapch', [], @isnumeric);
addOptional(p, 'rmvch', [], @isnumeric);
addOptional(p, 'fsIn', 20000, @isnumeric);
addOptional(p, 'accCh', [], @isnumeric);

parse(p, varargin{:})
basepath = p.Results.basepath;
exp = p.Results.exp;
rec = p.Results.rec;
nchans = p.Results.nchans;
mapch = p.Results.mapch;
rmvch = p.Results.rmvch;
fsIn = p.Results.fsIn;
accCh = p.Results.accCh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get files and folders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nprocessing %s\n\n', basepath)
% enter record node
cd(basepath)
recNodeDir = dir;
recNodeDir = recNodeDir(3);
recNodeName = recNodeDir.name;
cd(recNodeName);

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
    exPath{i} = fullfile(expFolders(expIdx(i)).folder, expNames{expIdx(i)});
    % get time from xml file
    if isequal(exp(i), 1)
        xmlname = ['settings.xml'];
    else
        xmlname = ['settings_' num2str(exp(i)) '.xml'];
    end
    t = convert(xmltree(xmlname));
    expTime{i} = datetime(t.INFO.DATE, 'InputFormat', 'dd MMM yyyy HH:mm:ss');
    
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
        recIdx = find(strcmp(recNames, ['recording' num2str(rec{exp(i)}(ii))]));
        recPath{k} = fullfile(recFolders(recIdx).folder, recNames{recIdx});
        
        % get recording start time
        recSyncFile = dir(recPath{k});
        recSyncFile = recSyncFile(contains({recSyncFile.name}, 'sync_messages'));
        recSyncFile = fullfile(recSyncFile.folder, recSyncFile.name);
        fid = fopen(recSyncFile, 'r');
        txt = fscanf(fid, '%s');
        txt = split(txt, ':');
        txt = split(txt{end}, '@');
        if ~isempty(txt{1})
            recStart(ii) = str2num(txt{1}) / fsIn;
        else
            recStart(ii) = 60 * 5;
        end
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
    % create new basefolder
    [newpath, baseTime] = fileparts(basepath);
    basename = bz_BasenameFromBasepath(newpath);
    dn = datenum(baseTime, 'yyyy-MM-dd');
    newpath = fullfile(newpath, [basename '_' datestr(dn, 'yyMMdd')]);
    expName = datestr(expTime{1} + seconds(recStart(1)), 'HHMMss');
    exPathNew = [newpath, '_', expName];
    [~, expName] = fileparts(exPathNew);
    mkdir(exPathNew)
    fprintf('created %s\n', exPathNew)
       
    % concatenate timestamps.npy and make sure dat files are not zero padded
    cat_OE_tstamps('orig_paths', recPath, 'new_path', exPathNew,...
        'nchans', nchans);
    
    % pre-process dat
    datInfo = preprocDat('orig_paths', recPath,...
        'newfile', fullfile(exPathNew, [expName, '.dat']), 'mapch', mapch,...
        'rmvch', rmvch, 'nchans', nchans, 'saveVar', true,...
        'chunksize', 5e6, 'precision', 'int16');  
    
    % get digital input
    getDinOE('basepath', recPath, 'newpath', exPathNew,...
        'concat', true, 'nchans', [], 'precision', 'int16',...
        'saveVar', true);
    
    % create lfp file
    LFPfromDat('basepath', exPathNew, 'cf', 450, 'chunksize', 5e6,...
        'nchans', length(mapch) - length(rmvch), 'fsOut', 1250,...
        'fsIn', fsIn)
    
    % get acceleration
    if ~isempty(accCh)
        EMGfromACC('basepath', exPathNew, 'fname', '',...
            'nchans', length(mapch) - length(rmvch), 'ch', accCh, 'force', false, 'saveVar', true,...
            'graphics', false, 'fsOut', 1250, 'fsIn', 1250);
    end
    
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