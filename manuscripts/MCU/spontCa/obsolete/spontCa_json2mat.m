function spontCa_json2mat(varargin)
% SPONTCA_JSON2MAT  Merge per-image LLM marks into per-cell event files.
%
% Reads <llmDir>/raw/<sbjID>_w*.json (one per rendered window, written
% by the image-marking subagents). Deduplicates events across
% overlapping windows, computes amp/dur/int from the trace, and writes
% one <llmDir>/<sbjID>.mat per cell as a bare events table with columns
% {compartment, start, stop, amp, dur, int}.
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
%   'outDir'   - dir to write <sbjID>.mat into. Default: <spontCa>/llm/.
%   'tolDedup' - cross-window dedup tolerance (s). Default 0.5.
%
% See also: SPONTCA_RENDER, SPONTCA_LLMCUR, SPONTCA_COMPARE,
%           SPONTCA_WRITEEVENTS

%% ARGUMENTS

p = inputParser;
p.addParameter('cells', {}, @(x) iscell(x) || ischar(x) || isstring(x));
p.addParameter('rawDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('outDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('tolDedup', 0.5, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0);
parse(p, varargin{:});
P = p.Results;
if ischar(P.cells) || isstring(P.cells), P.cells = cellstr(P.cells); end

thisDir = fileparts(mfilename('fullpath'));
llmDir = fullfile(thisDir, 'llm');
if isempty(P.rawDir), P.rawDir = fullfile(llmDir, 'raw'); end
if isempty(P.outDir), P.outDir = llmDir; end
P.rawDir = char(P.rawDir);
P.outDir = char(P.outDir);
if ~exist(P.outDir, 'dir'), mkdir(P.outDir); end


%% LOAD TRACES

[tblCell, fs] = spontCa_load();
dt = 1 / fs;
nT = size(tblCell.trace, 2);


%% CELL LIST FROM AVAILABLE RAW JSONS

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


%% ASSEMBLE PER CELL

for iCell = 1:numel(cellsToDo)
    sName = cellsToDo{iCell};
    iC = find(tblCell.sbjID == sName & tblCell.compartment == 'Cyto');
    iM = find(tblCell.sbjID == sName & tblCell.compartment == 'Mito');
    if isempty(iC) || isempty(iM)
        warning('Cell %s missing cyto or mito row; skipped', sName);
        continue;
    end
    cyTrace = tblCell.trace(iC, :);
    miTrace = tblCell.trace(iM, :);
    bslCyto = rollingPercentile(cyTrace, round(30 * fs), 20);
    bslMito = rollingPercentile(miTrace, round(30 * fs), 20);

    matches = startsWith(rawNames, [sName '_w']);
    windowFiles = rawNames(matches);
    [cyStarts, miStarts, miStops] = collectRawWindows( ...
        windowFiles, P.rawDir);

    [cyStarts, ~]       = dedupTimes(cyStarts, [], P.tolDedup);
    [miStarts, miStops] = dedupTimes(miStarts, miStops, P.tolDedup);

    cyEv = buildCytoEvents(cyStarts, cyTrace, bslCyto, fs, dt, nT);
    miEv = buildMitoEvents(miStarts, miStops, miTrace, bslMito, ...
        fs, dt, nT);
    events = [cyEv; miEv]; %#ok<NASGU>

    fpath = fullfile(P.outDir, sprintf('%s.mat', sName));
    if exist(fpath, 'file')
        bkupDir = fullfile(P.outDir, 'bkup');
        if ~exist(bkupDir, 'dir'), mkdir(bkupDir); end
        stamp = datestr(now, 'yymmdd_HHMMSS'); %#ok<TNOW1,DATST>
        copyfile(fpath, fullfile(bkupDir, ...
            sprintf('%s_%s.mat', sName, stamp)));
    end
    save(fpath, 'events');
    fprintf('[spontCa_json2mat] %s : %d cyto, %d mito -> %s\n', ...
        sName, height(cyEv), height(miEv), fpath);
end
end


%% HELPERS

function [cy, mi, mistops] = collectRawWindows(windowFiles, rawDir)
cy = []; mi = []; mistops = [];
for iF = 1:numel(windowFiles)
    f = fopen(fullfile(rawDir, windowFiles{iF}), 'r');
    txt = fread(f, Inf, '*char')';
    fclose(f);
    try
        mark = jsondecode(txt);
    catch ME
        warning('Failed to parse %s: %s', windowFiles{iF}, ME.message);
        continue;
    end
    if isfield(mark, 'cyto_peaks')
        cy = [cy; mark.cyto_peaks(:)]; %#ok<AGROW>
    end
    if isfield(mark, 'mito_peaks') && ~isempty(mark.mito_peaks)
        if isstruct(mark.mito_peaks)
            if numel(mark.mito_peaks) == 1 && ...
                    numel(mark.mito_peaks.start) > 1
                mi      = [mi;      mark.mito_peaks.start(:)]; %#ok<AGROW>
                mistops = [mistops; mark.mito_peaks.stop(:)];  %#ok<AGROW>
            else
                for k = 1:numel(mark.mito_peaks)
                    mi(end+1, 1)      = mark.mito_peaks(k).start; %#ok<AGROW>
                    mistops(end+1, 1) = mark.mito_peaks(k).stop;  %#ok<AGROW>
                end
            end
        elseif iscell(mark.mito_peaks)
            for k = 1:numel(mark.mito_peaks)
                mi(end+1, 1)      = mark.mito_peaks{k}.start; %#ok<AGROW>
                mistops(end+1, 1) = mark.mito_peaks{k}.stop;  %#ok<AGROW>
            end
        end
    end
end
end


function [s, e] = dedupTimes(s, e, tol)
if isempty(s), return; end
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


function tbl = buildCytoEvents(starts, trace, bsl, fs, dt, nT)
% Cyto stops are placeholders (start + 2 samples), dur = 2*dt, int = 0.
starts = starts(:);
starts = starts(starts >= 0 & starts <= (nT - 1) * dt);
n = numel(starts);
amp = zeros(n, 1);
for k = 1:n
    smp = max(1, min(nT, round(starts(k) * fs) + 1));
    amp(k) = trace(smp) - bsl(smp);
end
tbl = table( ...
    repmat(categorical({'Cyto'}, {'Cyto','Mito'}), n, 1), ...
    starts, starts + 2 * dt, amp, ...
    repmat(2 * dt, n, 1), zeros(n, 1), ...
    'VariableNames', {'compartment', 'start', 'stop', 'amp', 'dur', 'int'});
end


function tbl = buildMitoEvents(starts, stops, trace, bsl, fs, dt, nT)
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
tbl = table( ...
    repmat(categorical({'Mito'}, {'Cyto','Mito'}), n, 1), ...
    starts, stops, amp, stops - starts, intg, ...
    'VariableNames', {'compartment', 'start', 'stop', 'amp', 'dur', 'int'});
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
