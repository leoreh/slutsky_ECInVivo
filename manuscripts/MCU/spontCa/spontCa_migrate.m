function spontCa_migrate(varargin)
% SPONTCA_MIGRATE  One-shot migration from the old per-compartment
% curation files to the new single-file-per-cell events-table format.
%
% Reads <oldDir>/<sbjID>_Cyto.mat and <oldDir>/<sbjID>_Mito.mat (old
% format: struct cur with .start .stop .amp .dur .int) and writes
% <newDir>/<sbjID>.mat (new format: struct cur with .events table). Old
% files are NOT deleted; user can clean up spontCa_curated/ manually
% after verifying.
%
% USAGE
%   spontCa_migrate()                       % default paths
%   spontCa_migrate('oldDir', X, 'newDir', Y)
%   spontCa_migrate('overwrite', true)      % replace existing new-format
%
% OPTIONAL (Name-Value):
%   'oldDir'    - source dir (old format). Default: spontCa_curated/.
%   'newDir'    - target dir (new format). Default: man/.
%   'overwrite' - replace existing new-format files. Default false.

p = inputParser;
p.addParameter('oldDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('newDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('overwrite', false, @islogical);
parse(p, varargin{:});
P = p.Results;

thisDir = fileparts(mfilename('fullpath'));
if isempty(P.oldDir), P.oldDir = fullfile(thisDir, 'spontCa_curated'); end
if isempty(P.newDir), P.newDir = fullfile(thisDir, 'man'); end
P.oldDir = char(P.oldDir);
P.newDir = char(P.newDir);
if ~exist(P.oldDir, 'dir')
    error('Source dir not found: %s', P.oldDir);
end
if ~exist(P.newDir, 'dir'), mkdir(P.newDir); end

% Find <sbjID>_Cyto.mat / <sbjID>_Mito.mat pairs at the TOP level only
% (not bkup/). Pair them up by sbjID prefix.
files = dir(fullfile(P.oldDir, '*_Cyto.mat'));
nMigrated = 0;
for k = 1:numel(files)
    sName = regexprep(files(k).name, '_Cyto\.mat$', '');
    cyPath = fullfile(P.oldDir, [sName '_Cyto.mat']);
    miPath = fullfile(P.oldDir, [sName '_Mito.mat']);
    if ~exist(miPath, 'file')
        warning('No mito pair for %s; skipped', sName);
        continue;
    end
    outPath = fullfile(P.newDir, [sName '.mat']);
    if exist(outPath, 'file') && ~P.overwrite
        fprintf('[migrate] %s already in new format; skipped\n', sName);
        continue;
    end

    cy = load(cyPath, 'cur'); cyOld = cy.cur;
    mi = load(miPath, 'cur'); miOld = mi.cur;

    evCyto = struct( ...
        'start', cyOld.start(:), 'stop', cyOld.stop(:), ...
        'amp',   cyOld.amp(:),   'dur',  cyOld.dur(:), ...
        'int',   cyOld.int(:));
    evMito = struct( ...
        'start', miOld.start(:), 'stop', miOld.stop(:), ...
        'amp',   miOld.amp(:),   'dur',  miOld.dur(:), ...
        'int',   miOld.int(:));

    if isfield(cyOld, 'fs'), fs = cyOld.fs;
    elseif isfield(miOld, 'fs'), fs = miOld.fs;
    else, fs = NaN;
    end
    cur = struct( ...
        'sbjID',   sName, ...
        'fs',      fs, ...
        'savedAt', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ... %#ok<TNOW1,DATST>
        'source',  'man', ...
        'events',  [spontCa_ev2tbl(evCyto, 'Cyto'); ...
                    spontCa_ev2tbl(evMito, 'Mito')]);
    save(outPath, 'cur');
    nMigrated = nMigrated + 1;
    fprintf('[migrate] %s -> %s\n', sName, outPath);
end

fprintf('[migrate] %d cells migrated. Old files preserved in %s.\n', ...
    nMigrated, P.oldDir);
end
