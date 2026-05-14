function spontCa_writeEvents(tblEvent, dirPath, varargin)
% SPONTCA_WRITEEVENTS  Split a long-format tblEvent by sbjID and write
% one bare events table per cell into dirPath.
%
% Per-cell file: <sbjID>.mat containing a variable named `events`, a
% table with the minimal schema {compartment, start, stop}. No sbjID
% column (it's the filename), no metadata wrapper. Any of {amp, dur,
% int} in the input are dropped before write - these are recomputed on
% demand in spontCa2_metrics. For Cyto rows, stop is forced to equal
% start (events are single-sample at fs=3 Hz; baking this in keeps the
% on-disk schema unambiguous).
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

% Strip non-essential columns. amp/dur/int are recomputed downstream.
dropVars = intersect(tblEvent.Properties.VariableNames, ...
    {'amp', 'dur', 'int', 'flux'});
if ~isempty(dropVars)
    tblEvent = removevars(tblEvent, dropVars);
end

% Enforce single-sample cyto: stop == start.
if all(ismember({'compartment', 'start', 'stop'}, ...
        tblEvent.Properties.VariableNames))
    isCyto = tblEvent.compartment == 'Cyto';
    tblEvent.stop(isCyto) = tblEvent.start(isCyto);
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
