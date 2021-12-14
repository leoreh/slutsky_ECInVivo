function [sessionInfo] = getSessionInfo(basePath,varargin)
% [sessionInfo] = getSessionInfo(basePath) loads the sessionInfo metadata
% for the recording in basePath. basePath should be in the format:
%       /whateverPath/baseName/
%           a file  basePath/baseName.sessionInfo.mat
%           or      basePath/baseName.xml
%           should exist.
% If no baseName.sessionInfo.mat exists, loads from the xml.
%
% INPUT
%   basePath            directory: '/whatevetPath/baseName/'
%   (options)
%       'noPrompts'     (default: false) prevents prompts about
%                       saving/adding metadata
%       'editGUI'       (default: false) opens a GUI to edit select
%                       sessionInfo fields (beta, please add/improve!)
%
% OUTPUT
%   sessionInfo         metadata structure
%
% CALLS
%   LoadParameters      Originally from FMAToolbox, edited in buzlab
% 
% 2017 DLevenstein and DTingley
% 
% UPDATES
% 23 nov 18 LH - input basePath not include baseName
% 06 oct 19 LH - send baseName to LoadParameters
% 

%% inputs and defaults
p = inputParser;
addParameter(p,'noPrompts',false,@islogical);
addParameter(p,'editGUI',false,@islogical);
parse(p,varargin{:})
noPrompts = p.Results.noPrompts;
editGUI = p.Results.editGUI;

if ~exist('basePath','var')
    basePath = pwd;
end
if ispc
    [~, baseName] = fileparts(basePath);
else
    baseName = bz_BasenameFromBasepath(basePath);
end
filename = fullfile(basePath,[baseName,'.sessionInfo.mat']);

%% Load the stuff
%d = dir('*sessionInfo*'); %all files with sessioninfo in the name
if exist(filename,'file')
    sessionInfostruct = load(filename);
    %Checks that there is a single structure in the sessionInfo file
    varsInFile = fieldnames(sessionInfostruct); 
    if numel(varsInFile)==1
        sessionInfo = sessionInfostruct.(varsInFile{1});
    else
        warning('Your .sessionInfo.mat has multiple variables/structures in it... wtf.')
        sessionInfo = sessionInfostruct;
    end
    SIexist = true;  %Marks that session info exists as expected
else
   warning(['could not find file ',baseName,'.sessionInfo.mat ',...
       'running LoadParameters instead..']) 
   sessionInfo = LoadParameters([baseName '.xml']);
   SIexist = false; 
end

%% Check sessionInfo using bz_isSessionInfo and update if necessary
% bz_isSessionInfo(sessionInfo);

%Here: prompt user to add any missing sessionInfo fields and save
if editGUI
    [ sessionInfo ] = bz_sessionInfoGUI(sessionInfo);
    SIexist = false;
elseif ~isfield(sessionInfo,'region') && ~noPrompts
    regionadd = questdlg(['Your sessionInfo is missing regions, ',...
        'would you like to add them?'],'Add Regions?','Yes');
    switch regionadd
        case 'Yes'
            [sessionInfo] = bz_sessionInfoGUI(sessionInfo,'Regions');
            SIexist = false; 
        case 'Cancel'
            return
    end
end

%Should check that sessionInfo.session.name and sesioninfo.session.path
%match basePath....  if not, prompt the user to save with the correct
%information.
    
%% Save sessionInfo file   
%Prompt user to save basePath/baseName.sessionInfo.mat 
%if loaded from xml or changed
if ~noPrompts && ~SIexist %Inform the user that they should save a file for later
    savebutton = questdlg(['Would you like to save your sessionInfo in ',...
        filename '?'],'Save sessionInfo?','Yes');
    switch savebutton
        case 'Yes'
            save(filename,'sessionInfo'); 
        case 'Cancel'
            return
    end
end

end
