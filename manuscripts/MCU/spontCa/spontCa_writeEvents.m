function spontCa_writeEvents(tblEvent, dirPath, varargin)
% SPONTCA_WRITEEVENTS  Split a long-format tblEvent by sbjID and write
% one bare events table per cell into dirPath.
%
% Per-cell file: <sbjID>.mat containing a variable named `events`, a
% table with columns {compartment, start, stop, amp, dur, int}. No
% sbjID column (it's the filename), no metadata wrapper.
%
% If a target file exists and 'backup' is true (default), the existing
% file is copied to <dirPath>/bkup/<sbjID>_<stamp>.mat before overwrite.
%
% USAGE
%   spontCa_writeEvents(tblEvent, autoDir)
%   spontCa_writeEvents(tblEvent, manDir, 'backup', false)
%
% OPTIONAL (Name-Value):
%   'backup' - back up existing files to dirPath/bkup/. Default true.
%
% See also: SPONTCA_READEVENTS, SPONTCA_DETECT

p = inputParser;
addRequired(p, 'tblEvent', @istable);
addRequired(p, 'dirPath',  @(x) ischar(x) || isstring(x));
addParameter(p, 'backup', true, @islogical);
parse(p, tblEvent, dirPath, varargin{:});
P = p.Results;

dirPath = char(P.dirPath);
if ~exist(dirPath, 'dir'), mkdir(dirPath); end

if ~ismember('sbjID', tblEvent.Properties.VariableNames)
    error('spontCa_writeEvents:noSbjID', ...
        'tblEvent must have an sbjID column');
end

sbjIDs = unique(tblEvent.sbjID, 'stable');
stamp = datestr(now, 'yymmdd_HHMMSS'); %#ok<TNOW1,DATST>

for k = 1:numel(sbjIDs)
    sid = sbjIDs(k);
    sub = tblEvent(tblEvent.sbjID == sid, :);
    sub = removevars(sub, 'sbjID');
    events = sub; %#ok<NASGU>
    fpath = fullfile(dirPath, sprintf('%s.mat', char(sid)));
    if P.backup && exist(fpath, 'file')
        bkupDir = fullfile(dirPath, 'bkup');
        if ~exist(bkupDir, 'dir'), mkdir(bkupDir); end
        copyfile(fpath, fullfile(bkupDir, ...
            sprintf('%s_%s.mat', char(sid), stamp)));
    end
    save(fpath, 'events');
end
end
