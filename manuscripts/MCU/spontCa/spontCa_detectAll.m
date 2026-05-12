function tbl = spontCa_detectAll(tbl, fs, varargin)
% SPONTCA_DETECTALL  Run spontCa_detect on every row of a long-format
% tbl and write one auto/<sbjID>.mat per cell.
%
% Replaces the inline detection loop in mcu_spontCa. The file format
% matches what spontCa_manCur and spontCa_json2mat write: a struct cur
% with {sbjID, fs, savedAt, source='auto', events table}. The same
% three folders (auto/, man/, llm/) hold interchangeable per-cell
% event files.
%
% Adds per-row event columns to tbl (start/stop/amp/dur/int as cell
% arrays of column vectors) so downstream code (spontCa_manCur,
% spontCa_finalize) can read them without round-tripping disk.
%
% USAGE
%   tbl = spontCa_detectAll(tbl, fs);
%   tbl = spontCa_detectAll(tbl, fs, ...
%               'paramsCyto', {'minAmp', 0.05, 'minIEI', 1.0, ...}, ...
%               'paramsMito', {...});
%
% OPTIONAL (Name-Value):
%   'paramsCyto' - cell array of NV pairs forwarded to spontCa_detect
%                  for Cyto rows. Default {}.
%   'paramsMito' - cell array of NV pairs forwarded to spontCa_detect
%                  for Mito rows. Default {}.
%   'autoDir'    - output directory. Default <spontCa>/auto/.
%   'overwrite'  - re-write existing auto/<sbjID>.mat files. Default
%                  true (autodetection is cheap; stale auto files would
%                  bias downstream comparisons).
%   'verbose'    - print per-cell summary. Default true.
%
% See also: SPONTCA_DETECT, SPONTCA_EV2TBL, SPONTCA_MANCUR,
%           SPONTCA_FINALIZE

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'fs',  @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'paramsCyto', {}, @iscell);
addParameter(p, 'paramsMito', {}, @iscell);
addParameter(p, 'autoDir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'overwrite', true, @islogical);
addParameter(p, 'verbose', true, @islogical);
parse(p, tbl, fs, varargin{:});
P = p.Results;

if isempty(P.autoDir)
    P.autoDir = fullfile(fileparts(mfilename('fullpath')), 'auto');
end
P.autoDir = char(P.autoDir);
if ~exist(P.autoDir, 'dir'), mkdir(P.autoDir); end

n = height(tbl);
tbl.start = cell(n, 1);
tbl.stop  = cell(n, 1);
tbl.amp   = cell(n, 1);
tbl.dur   = cell(n, 1);
tbl.int   = cell(n, 1);

% Per-trace detection
for iRow = 1:n
    if tbl.compartment(iRow) == 'Cyto'
        ev = spontCa_detect(tbl.trace(iRow, :), fs, P.paramsCyto{:});
    else
        ev = spontCa_detect(tbl.trace(iRow, :), fs, P.paramsMito{:});
    end
    tbl.start{iRow} = ev.start(:);
    tbl.stop{iRow}  = ev.stop(:);
    tbl.amp{iRow}   = ev.amp(:);
    tbl.dur{iRow}   = ev.dur(:);
    tbl.int{iRow}   = ev.int(:);
end

% Per-cell file write
cells = unique(cellstr(string(tbl.sbjID)), 'stable');
for iCell = 1:numel(cells)
    sid = cells{iCell};
    outPath = fullfile(P.autoDir, [sid '.mat']);
    if exist(outPath, 'file') && ~P.overwrite, continue; end

    iC = find(tbl.sbjID == sid & tbl.compartment == 'Cyto');
    iM = find(tbl.sbjID == sid & tbl.compartment == 'Mito');
    evCyto = struct('start', tbl.start{iC}, 'stop', tbl.stop{iC}, ...
        'amp', tbl.amp{iC}, 'dur', tbl.dur{iC}, 'int', tbl.int{iC});
    evMito = struct('start', tbl.start{iM}, 'stop', tbl.stop{iM}, ...
        'amp', tbl.amp{iM}, 'dur', tbl.dur{iM}, 'int', tbl.int{iM});
    cur = struct( ...
        'sbjID',   sid, ...
        'fs',      fs, ...
        'savedAt', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ... %#ok<TNOW1,DATST>
        'source',  'auto', ...
        'events',  [spontCa_ev2tbl(evCyto, 'Cyto'); ...
                    spontCa_ev2tbl(evMito, 'Mito')]);
    save(outPath, 'cur');
end

if P.verbose
    nCyto = sum(cellfun(@numel, tbl.start(tbl.compartment == 'Cyto')));
    nMito = sum(cellfun(@numel, tbl.start(tbl.compartment == 'Mito')));
    fprintf('[spontCa_detectAll] %d cells | %d cyto + %d mito events | %s\n', ...
        numel(cells), nCyto, nMito, P.autoDir);
end

end
