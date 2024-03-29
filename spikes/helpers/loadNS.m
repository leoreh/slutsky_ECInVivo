function [data, ndata] = loadNS(varargin)

% loads data stored in neurosuite format (clu, res, fet, or spk).
% pwd must be path to folder. fs, nfet, and spkgrp can be extracted from
% session info if exists
%
% INPUT:
%   basepath    string. path to recording folder {pwd}.
%   basename    string. basename of files. if empty will be extracted from
%               basepath
%   datatype    string. data to load; 'res', 'spk', {'clu'}, or 'fet'
%   grpid        numeric. spkgrp number {1}
%   fs          numeric. sample frequency. if empty will try to take from
%               session (cell explorer format)
%   nfet        numeric. number of features to keep from fet data {12}
%   spkgrp      array where each cell is the electrodes for a spike group
%   session     struct. session info (see CE_sessionTemplate)
%   nspks       numeric. given for spk. if empty will extract from clu. 
% 
% OUTPUT:
%   data        clu / res / spk / fet
%   ndata       nclu (if clu) or nfet (if fet)
%
% DEPENDENCIES
%
% TO DO LIST
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
addOptional(p, 'fs', [], @isnumeric);
addOptional(p, 'nfet', [], @isnumeric);
addOptional(p, 'spkgrp', {}, @iscell);
addOptional(p, 'session', []);
addOptional(p, 'nspks', [], @isnumeric);

parse(p, varargin{:})
basepath    = p.Results.basepath;
basename    = p.Results.basename;
datatype    = p.Results.datatype;
grpid       = p.Results.grpid;
fs          = p.Results.fs;
nfet        = p.Results.nfet;
spkgrp      = p.Results.spkgrp;
session     = p.Results.session;
nspks       = p.Results.nspks;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(session)
    sessionName = [basename, '.session.mat'];
    if ~exist(sessionName, 'file')
        session = CE_sessionTemplate(pwd, 'viaGUI', false,...
            'forceL', true, 'saveVar', true);
    else
        load(sessionName)
    end
end
if isempty(basename)
    [~, basename] = fileparts(basepath);
end
if isempty(spkgrp)
    spkgrp = session.extracellular.spikeGroups.channels;
end
filename = sprintf('%s.%s.%d', basename, datatype, grpid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('loading %s... ', filename)
switch datatype
    case 'clu'    
        fid = fopen(filename, 'r');
        ndata = fscanf(fid, '%d\n', 1);
        data = fscanf(fid, '%d\n');
        rc = fclose(fid);
        uclu = unique(data);
        
    case 'res'
        fid = fopen(filename, 'r');
        data = fscanf(fid, '%d\n');
        rc = fclose(fid);
        
    case 'fet'      % note only features are returned (not energy etc.)
        if isempty(nfet)
            nfet = 3 * length(spkgrp{grpid});
        end
        fid = fopen(filename, 'r');
        ndata = fscanf(fid, '%d', 1);
        data = fscanf(fid, '%d', [ndata, inf])';
        data = data(:, 1 : nfet);
        rc = fclose(fid);
        
    case 'spk'
        if isempty(nspks)
            clu = loadNS('datatype', 'clu', 'grpid', grpid);
            nspks = length(clu);
        end
        if isempty(fs)
            fs = session.extracellular.sr;
        end
        sniplength = ceil(1.6 * 10^-3 * fs);
        fid = fopen(filename, 'r');
        data = fread(fid, 'int16');
        data = reshape(data, length(spkgrp{grpid}), sniplength, nspks);
        rc = fclose(fid);
end
if rc ~= 0 || isempty(data)
    warning('failed to read %s.%d\n',...
        datatype, num2str(grpid))
end
fprintf('done in %.2f s\n', toc)
end
% EOF