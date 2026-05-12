function spontCa_json2mat(varargin)
% SPONTCA_JSON2MAT  Merge per-image LLM marks into per-cell event files.
%
% Reads <llmDir>/raw/<sbjID>_w*.json (one per rendered window, written
% by the image-marking subagents). Deduplicates events across
% overlapping windows, computes amp/dur/int from the trace, and writes
% one <llmDir>/<sbjID>.mat per cell.
%
% File format (matches spontCa_manCur save format):
%   cur.sbjID    char
%   cur.fs       double
%   cur.savedAt  char timestamp
%   cur.source   'llm'
%   cur.events   table with columns
%                {compartment, start, stop, amp, dur, int}
%                - one row per event, sorted by (compartment, start)
%                - compartment is categorical {Cyto, Mito}
%
% Cyto stops are placeholders (start + 2 samples; cyto stops are not
% biologically real at fs=3); mito stops come from the LLM marks.
%
% USAGE
%   spontCa_json2mat()                       % all cells with raw JSONs
%   spontCa_json2mat('cells', {'Ctrl_03'})   % subset
%   spontCa_json2mat('outDir', dirPath)      % override output dir
%
% OPTIONAL (Name-Value):
%   'cells'    - cellstr of sbjIDs. Default: all with raw/*.json.
%   'rawDir'   - dir holding subagent JSONs.
%                Default: <llmDir>/raw/ where llmDir = <spontCa>/llm/.
%   'outDir'   - dir to write <sbjID>.mat into.
%                Default: <spontCa>/llm/.
%   'tolDedup' - cross-window dedup tolerance (s). Default 0.5.
%
% See also: SPONTCA_RENDER, SPONTCA_LLMCUR, SPONTCA_COMPARE,
%           SPONTCA_FINALIZE, SPONTCA_MANCUR, SPONTCA_EV2TBL

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
p.addParameter('cells', {}, @(x) iscell(x) || ischar(x) || isstring(x));
p.addParameter('rawDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('outDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('tolDedup', 0.5, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0);
parse(p, varargin{:});
P = p.Results;
if ischar(P.cells) || isstring(P.cells)
    P.cells = cellstr(P.cells);
end

thisDir = fileparts(mfilename('fullpath'));
llmDir = fullfile(thisDir, 'llm');
if isempty(P.rawDir), P.rawDir = fullfile(llmDir, 'raw'); end
if isempty(P.outDir), P.outDir = llmDir; end
P.rawDir = char(P.rawDir);
P.outDir = char(P.outDir);
if ~exist(P.outDir, 'dir'), mkdir(P.outDir); end


%% ========================================================================
%  LOAD TRACES
%  ========================================================================

[tbl, fs] = spontCa_load();
dt = 1 / fs;
nT = size(tbl.trace, 2);


%% ========================================================================
%  CELL LIST FROM AVAILABLE RAW JSONS
%  ========================================================================

allRaw = dir(fullfile(P.rawDir, '*_w*.json'));
if isempty(allRaw)
    error('No raw JSONs in %s', P.rawDir);
end
rawNames = {allRaw.name};
sbjFromName = regexprep(rawNames, '_w\d+\.json$', '');
rawCells    = unique(sbjFromName, 'stable');

if isempty(P.cells)
    cellsToDo = rawCells;
else
    cellsToDo = intersect(rawCells, P.cells, 'stable');
end


%% ========================================================================
%  ASSEMBLE PER CELL
%  ========================================================================

for iCell = 1:numel(cellsToDo)
    sName = cellsToDo{iCell};
    iC = find(tbl.sbjID == sName & tbl.compartment == 'Cyto');
    iM = find(tbl.sbjID == sName & tbl.compartment == 'Mito');
    if isempty(iC) || isempty(iM)
        warning('Cell %s missing cyto or mito row; skipped', sName);
        continue;
    end
    cyTrace = tbl.trace(iC, :);
    miTrace = tbl.trace(iM, :);
    bslCyto = rollingPercentile(cyTrace, round(30 * fs), 20);
    bslMito = rollingPercentile(miTrace, round(30 * fs), 20);

    matches = startsWith(rawNames, [sName '_w']);
    windowFiles = rawNames(matches);

    cyStarts = [];
    miStarts = [];
    miStops  = [];

    for iF = 1:numel(windowFiles)
        f = fopen(fullfile(P.rawDir, windowFiles{iF}), 'r');
        txt = fread(f, Inf, '*char')';
        fclose(f);
        try
            mark = jsondecode(txt);
        catch ME
            warning('Failed to parse %s: %s', windowFiles{iF}, ME.message);
            continue;
        end
        if isfield(mark, 'cyto_peaks')
            cyStarts = [cyStarts; mark.cyto_peaks(:)]; %#ok<AGROW>
        end
        if isfield(mark, 'mito_peaks') && ~isempty(mark.mito_peaks)
            if isstruct(mark.mito_peaks)
                if numel(mark.mito_peaks) == 1 && ...
                        numel(mark.mito_peaks.start) > 1
                    miStarts = [miStarts; mark.mito_peaks.start(:)]; %#ok<AGROW>
                    miStops  = [miStops;  mark.mito_peaks.stop(:)];  %#ok<AGROW>
                else
                    for k = 1:numel(mark.mito_peaks)
                        miStarts(end+1, 1) = mark.mito_peaks(k).start; %#ok<AGROW>
                        miStops(end+1, 1)  = mark.mito_peaks(k).stop;  %#ok<AGROW>
                    end
                end
            elseif iscell(mark.mito_peaks)
                for k = 1:numel(mark.mito_peaks)
                    miStarts(end+1, 1) = mark.mito_peaks{k}.start; %#ok<AGROW>
                    miStops(end+1, 1)  = mark.mito_peaks{k}.stop;  %#ok<AGROW>
                end
            end
        end
    end

    [cyStarts, ~]       = dedupTimes(cyStarts, [], P.tolDedup);
    [miStarts, miStops] = dedupTimes(miStarts, miStops, P.tolDedup);

    evCyto = buildCytoEv(cyStarts, cyTrace, bslCyto, fs, dt, nT);
    evMito = buildMitoEv(miStarts, miStops, miTrace, bslMito, fs, dt, nT);

    events = [spontCa_ev2tbl(evCyto, 'Cyto'); ...
              spontCa_ev2tbl(evMito, 'Mito')];

    cur = struct( ...
        'sbjID',   sName, ...
        'fs',      fs, ...
        'savedAt', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ... %#ok<TNOW1,DATST>
        'source',  'llm', ...
        'events',  events);

    outPath = fullfile(P.outDir, sprintf('%s.mat', sName));
    bkupExisting(outPath, P.outDir);
    save(outPath, 'cur');

    fprintf('[spontCa_json2mat] %s : %d cyto, %d mito -> %s\n', ...
        sName, height(events(events.compartment == 'Cyto', :)), ...
        height(events(events.compartment == 'Mito', :)), outPath);
end

end     % EOF


%% ========================================================================
%  HELPERS
%  ========================================================================

function [s, e] = dedupTimes(s, e, tol)
if isempty(s)
    return;
end
[s, ord] = sort(s);
if ~isempty(e), e = e(ord); end
keep = true(numel(s), 1);
for k = 2:numel(s)
    if s(k) - s(find(keep(1:k-1), 1, 'last')) < tol
        keep(k) = false;
    end
end
s = s(keep);
if ~isempty(e), e = e(keep); end
end


function ev = buildCytoEv(starts, trace, bsl, fs, dt, nT)
% Cyto: stop = start + 2 samples (placeholder), dur = 2*dt, int = 0.
starts = starts(:);
starts = starts(starts >= 0 & starts <= (nT - 1) * dt);
n = numel(starts);
amp = zeros(n, 1);
for k = 1:n
    smp = max(1, min(nT, round(starts(k) * fs) + 1));
    amp(k) = trace(smp) - bsl(smp);
end
ev = struct( ...
    'start', starts, ...
    'stop',  starts + 2 * dt, ...
    'amp',   amp, ...
    'dur',   repmat(2 * dt, n, 1), ...
    'int',   zeros(n, 1));
end


function ev = buildMitoEv(starts, stops, trace, bsl, fs, dt, nT)
starts = starts(:);
stops  = stops(:);
valid = starts >= 0 & starts <= (nT - 1) * dt & ...
        stops  >= starts & stops <= (nT - 1) * dt;
starts = starts(valid);
stops  = stops(valid);
n = numel(starts);
amp = zeros(n, 1);
intg = zeros(n, 1);
for k = 1:n
    pkSmp = max(1, min(nT, round(starts(k) * fs) + 1));
    stSmp = max(pkSmp, min(nT, round(stops(k) * fs) + 1));
    amp(k) = trace(pkSmp) - bsl(pkSmp);
    seg = trace(pkSmp:stSmp) - bsl(pkSmp:stSmp);
    seg(isnan(seg)) = 0;
    intg(k) = trapz(seg) * dt;
end
ev = struct( ...
    'start', starts, ...
    'stop',  stops, ...
    'amp',   amp, ...
    'dur',   stops - starts, ...
    'int',   intg);
end


function bkupExisting(fpath, outDir)
if ~exist(fpath, 'file'), return; end
bkupDir = fullfile(outDir, 'bkup');
if ~exist(bkupDir, 'dir'), mkdir(bkupDir); end
[~, base, ext] = fileparts(fpath);
stamp = datestr(now, 'yymmdd_HHMMSS'); %#ok<TNOW1,DATST>
copyfile(fpath, fullfile(bkupDir, sprintf('%s_%s%s', base, stamp, ext)));
end


function out = rollingPercentile(x, winSmp, q)
nT = length(x);
out = nan(1, nT);
halfWin = floor(winSmp / 2);
for i = 1:nT
    lo = max(1, i - halfWin);
    hi = min(nT, i + halfWin);
    seg = x(lo:hi);
    seg = seg(~isnan(seg));
    if ~isempty(seg)
        out(i) = prctile(seg, q);
    end
end
end
