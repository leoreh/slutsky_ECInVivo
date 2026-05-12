function llmCur_assemble(varargin)
% LLMCUR_ASSEMBLE  Merge per-image LLM marks into per-cell curated .mat
% files, in the same format spontCa_manCur writes.
%
% Reads llmCur/raw/<sbjID>_w*.json (one per rendered window, written by
% the image-marking subagents). Deduplicates events across overlapping
% windows, computes amp/dur/int from the trace, and writes
% spontCa_curated/<sbjID>_<cmp>.mat. Cyto stops are forced to
% start + 2 samples (placeholder; cyto stops are not biologically real
% at fs=3); mito stops come straight from the LLM marks.
%
% USAGE
%   llmCur_assemble()                       % all cells with raw JSONs
%   llmCur_assemble('cells', {'Ctrl_03'})   % subset
%   llmCur_assemble('outDir', dirPath)      % override curated output dir
%
% OPTIONAL (Name-Value):
%   'cells'      - cellstr of sbjIDs to assemble. Default: all with raw.
%   'rawDir'     - dir holding subagent JSONs.
%                  Default: <this_file_dir>/raw/.
%   'outDir'     - dir to write <sbjID>_<cmp>.mat into.
%                  Default: ../spontCa_curated/ relative to this file.
%   'tolDedup'   - cross-window dedup tolerance (s). Default 0.5.
%
% See also: LLMCUR_RENDER, LLMCUR_RUN, SPONTCA_FINALIZE, SPONTCA_MANCUR

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
if isempty(P.rawDir)
    P.rawDir = fullfile(thisDir, 'raw');
end
if isempty(P.outDir)
    P.outDir = fullfile(fileparts(thisDir), 'spontCa_curated');
end
P.rawDir = char(P.rawDir);
P.outDir = char(P.outDir);
if ~exist(P.outDir, 'dir'), mkdir(P.outDir); end


%% ========================================================================
%  LOAD TRACES
%  ========================================================================

spontCaDir = fileparts(thisDir);
if exist(spontCaDir, 'dir') && ~contains(lower(path), lower(spontCaDir))
    addpath(spontCaDir);
end
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
% Extract sbjID = substring before "_w<NN>.json"
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

    % Gather raw JSONs for this cell
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
            % Two encodings: array of structs {start, stop} or struct
            % with fields start/stop as vectors.
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

    % Dedup
    [cyStarts, ~] = dedupTimes(cyStarts, [], P.tolDedup);
    [miStarts, miStops] = dedupTimes(miStarts, miStops, P.tolDedup);

    % Build curated structs
    curCyto = buildCyto(sName, cyStarts, cyTrace, bslCyto, fs, dt, nT);
    curMito = buildMito(sName, miStarts, miStops, miTrace, bslMito, ...
        fs, dt, nT);

    % Save
    cur = curCyto;
    save(fullfile(P.outDir, sprintf('%s_Cyto.mat', sName)), 'cur');
    cur = curMito;
    save(fullfile(P.outDir, sprintf('%s_Mito.mat', sName)), 'cur');

    fprintf('[llmCur_assemble] %s : %d cyto, %d mito events -> %s\n', ...
        sName, numel(curCyto.start), numel(curMito.start), P.outDir);
end

end     % EOF


%% ========================================================================
%  HELPERS
%  ========================================================================

function [s, e] = dedupTimes(s, e, tol)
% Sort by start time; merge entries within tol seconds (keep first).
if isempty(s)
    return;
end
[s, ord] = sort(s);
if ~isempty(e)
    e = e(ord);
end
keep = true(numel(s), 1);
for k = 2:numel(s)
    if s(k) - s(find(keep(1:k-1), 1, 'last')) < tol
        keep(k) = false;
    end
end
s = s(keep);
if ~isempty(e), e = e(keep); end
end


function cur = buildCyto(sName, starts, trace, bsl, fs, dt, nT)
% Cyto: stop is placeholder (start + 2 samples), dur/int zeroed.
starts = starts(:);
starts = starts(starts >= 0 & starts <= (nT - 1) * dt);
n = numel(starts);
stop = starts + 2 * dt;
amp = zeros(n, 1);
for k = 1:n
    smp = max(1, min(nT, round(starts(k) * fs) + 1));
    amp(k) = trace(smp) - bsl(smp);
end
cur = struct( ...
    'sbjID',       sName, ...
    'compartment', 'Cyto', ...
    'fs',          fs, ...
    'start',       starts, ...
    'stop',        stop, ...
    'amp',         amp, ...
    'dur',         repmat(2 * dt, n, 1), ...
    'int',         zeros(n, 1), ...
    'savedAt',     datestr(now, 'yyyy-mm-dd HH:MM:SS')); %#ok<TNOW1,DATST>
end


function cur = buildMito(sName, starts, stops, trace, bsl, fs, dt, nT)
% Mito: real start AND stop, compute amp/dur/int.
starts = starts(:);
stops  = stops(:);
valid = starts >= 0 & starts <= (nT - 1) * dt & ...
        stops  >= starts & stops <= (nT - 1) * dt;
starts = starts(valid);
stops  = stops(valid);
n = numel(starts);
amp = zeros(n, 1);
intg = zeros(n, 1);
dur  = stops - starts;
for k = 1:n
    pkSmp = max(1, min(nT, round(starts(k) * fs) + 1));
    stSmp = max(pkSmp, min(nT, round(stops(k) * fs) + 1));
    amp(k) = trace(pkSmp) - bsl(pkSmp);
    seg = trace(pkSmp:stSmp) - bsl(pkSmp:stSmp);
    seg(isnan(seg)) = 0;
    intg(k) = trapz(seg) * dt;
end
cur = struct( ...
    'sbjID',       sName, ...
    'compartment', 'Mito', ...
    'fs',          fs, ...
    'start',       starts, ...
    'stop',        stops, ...
    'amp',         amp, ...
    'dur',         dur, ...
    'int',         intg, ...
    'savedAt',     datestr(now, 'yyyy-mm-dd HH:MM:SS')); %#ok<TNOW1,DATST>
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
