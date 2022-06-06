function saveNS(data, varargin)

% currently only clu is ready
%
% INPUT:
%   basepath    string. path to recording folder {pwd}.
%   basename    string. basename of files. if empty will be extracted from
%               basepath
%   data        numeric. data to save (clu / res / spk / fet)
%   datatype    string. data to load; 'res', 'spk', {'clu'}, or 'fet'
%   grpid       numeric. spkgrp number {1}
%   session     struct. session info (see CE_sessionTemplate)
%   nfet        numeric. number of columns in fet file {17}
%   spkgrp      array where each cell is the electrodes for a spike group
%   nspks       numeric. given for spk. if empty will extract from clu.
%   bkup        logical. create backup of file or not {true}
%   bkpath      string. path for backup files
%
% DEPENDENCIES
%
% TO DO LIST
%   # add option to calc PCA from spk for fet
%
% 22 mar 21 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'basename', '', @ischar);
addOptional(p, 'datatype', '', @ischar);
addOptional(p, 'grpid', 1, @isnumeric);
addOptional(p, 'spkgrp', {}, @iscell);
addOptional(p, 'session', []);
addOptional(p, 'bkup', true, @islogical);
addOptional(p, 'bkpath', '', @ischar);

parse(p, varargin{:})
basepath    = p.Results.basepath;
basename    = p.Results.basename;
datatype    = p.Results.datatype;
grpid        = p.Results.grpid;
spkgrp      = p.Results.spkgrp;
session     = p.Results.session;
bkup        = p.Results.bkup;
bkpath      = p.Results.bkpath;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
if isempty(basename)
    [~, basename] = fileparts(basepath);
end
if isempty(session)
    sessionName = [basename, '.session.mat'];
    if ~exist(sessionName, 'file')
        % session params
        session = CE_sessionTemplate(pwd, 'viaGUI', false,...
            'forceDef', true, 'forceL', false, 'saveVar', false);
    else
        load(sessionName)
    end
end
if isempty(spkgrp)
    spkgrp = session.extracellular.spikeGroups.channels;
end
filename = sprintf('%s.%s.%d', basename, datatype, grpid);
fprintf('saving %s... ', filename)

% bkup
if bkup
    source = fullfile(basepath, filename);
    if exist(source, 'file')
        if isempty(bkpath)
            bkpath = fullfile(basepath, 'kk', 'bkup');
        end
        mkdir(bkpath)
        destination = fullfile(bkpath, [filename, '.' datestr(datetime, 'ddmmyy_HHMMss')]);
        fprintf('creating backup... ')
        copyfile(source, destination)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch datatype
    case 'clu'
        fid = fopen(filename, 'w');
        fprintf(fid, '%d\n', length(unique(data)));
        fprintf(fid, '%d\n', data);
        rc = fclose(fid);
        
    case 'res'
        fid = fopen(filename, 'w');
        fprintf(fid, '%d\n', data);
        rc = fclose(fid);
        
    case 'fet'   
        nfet = size(data, 1);
        fid = fopen(filename, 'w');
        formatstring = '%d';
        for ifet = 2 : nfet
            formatstring = [formatstring, '\t%d'];
        end
        formatstring = [formatstring, '\n'];
        fprintf(fid, '%d\n', nfet);
        fprintf(fid, formatstring, data);
        rc = fclose(fid);
    
    case 'spk'
        fid = fopen(filename, 'w');
        fwrite(fid, data(:), 'int16');
        rc = fclose(fid);
end

if rc ~= 0 || isempty(data)
    warning('failed to write %s.%d\n',...
        datatype, num2str(grpid))
end
fprintf('done in %.2f s\n', toc)
end
% EOF