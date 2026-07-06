function bkPath = backup_file(file, varargin)
% BACKUP_FILE  Copy a file to a timestamped backup before it is overwritten.
%
%   bkPath = BACKUP_FILE(file) copies FILE to <dir>/bkup/<name>_<stamp><ext>,
%   where <dir> is the file's folder, <stamp> is datestr(now, 'yymmdd_HHMMSS'),
%   and <name>/<ext> are the file's name and extension. The bkup subfolder is
%   created if missing. If FILE does not exist yet, nothing is copied and bkPath
%   is '' (a first save has no original to preserve).
%
%   This is the shared form of the backup convention used inline in
%   spontCa_writeEvents, utypes_push and saveNS: keep every overwritten original
%   as a dated copy so a save can always be rolled back.
%
%   INPUTS:
%       file     - (Char) Path to the file about to be overwritten.
%       varargin - Parameter/Value pairs:
%           'subdir' - (Char) Backup subfolder name. {'bkup'}
%
%   OUTPUT:
%       bkPath   - (Char) Path of the backup copy, or '' if FILE was absent.
%
%   HISTORY:
%       05 Jul 2026 - factored from spontCa_writeEvents for guiPath_curate saves.

p = inputParser;
addRequired(p, 'file', @(x) ischar(x) || isstring(x));
addParameter(p, 'subdir', 'bkup', @(x) ischar(x) || isstring(x));
parse(p, file, varargin{:});
file   = char(p.Results.file);
subdir = char(p.Results.subdir);

bkPath = '';
if ~isfile(file), return; end

[dirPath, name, ext] = fileparts(file);
bkupDir = fullfile(dirPath, subdir);
if ~exist(bkupDir, 'dir'), mkdir(bkupDir); end
stamp  = datestr(now, 'yymmdd_HHMMSS'); %#ok<TNOW1,DATST>
bkPath = fullfile(bkupDir, sprintf('%s_%s%s', name, stamp, ext));
copyfile(file, bkPath);
end
