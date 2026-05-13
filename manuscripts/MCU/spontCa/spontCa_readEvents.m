function tblEvent = spontCa_readEvents(dirPath)
% SPONTCA_READEVENTS  Concatenate all per-cell event files in dirPath
% into a single long-format tblEvent.
%
% Each <sbjID>.mat in dirPath holds a bare events table (variable name
% `events`) with columns {compartment, start, stop, amp, dur, int}.
% This function tags each row with its sbjID (from the filename) and
% vertcats into one tblEvent with columns
% {sbjID, compartment, start, stop, amp, dur, int}.
%
% Subdirectories (bkup/, images/, raw/) are ignored.
%
% USAGE
%   tblEvent = spontCa_readEvents(manDir)
%   tblEvent = spontCa_readEvents(llmDir)
%
% Returns an empty (but correctly-typed) tblEvent if dirPath does not
% exist or has no .mat files.
%
% See also: SPONTCA_WRITEEVENTS, SPONTCA_DETECT, SPONTCA_FINALIZE

dirPath = char(dirPath);
files = dir(fullfile(dirPath, '*.mat'));
if isempty(files)
    tblEvent = emptyTblEvent();
    return;
end

canonVars = {'sbjID', 'compartment', 'start', 'stop', 'amp', 'dur', 'int'};
chunks = cell(numel(files), 1);
for k = 1:numel(files)
    [~, sid] = fileparts(files(k).name);
    fpath = fullfile(dirPath, files(k).name);
    S = load(fpath);
    if isfield(S, 'events')
        events = S.events;
    elseif isfield(S, 'cur') && isfield(S.cur, 'events')
        % Legacy v4 format: struct cur wrapping the events table.
        events = S.cur.events;
    else
        warning('spontCa_readEvents:noTable', ...
            'Skipped %s: no events variable', files(k).name);
        continue;
    end
    nE = height(events);
    if nE > 0
        events.sbjID = repmat(categorical({sid}), nE, 1);
    else
        events.sbjID = categorical(strings(0, 1));
    end
    % Canonicalize: drop unexpected vars, fill missing ones with NaN/<undefined>
    missing = setdiff(canonVars, events.Properties.VariableNames);
    for v = missing
        if any(strcmp(v{1}, {'sbjID', 'compartment'}))
            events.(v{1}) = categorical(strings(nE, 1));
        else
            events.(v{1}) = nan(nE, 1);
        end
    end
    events = events(:, canonVars);
    chunks{k} = events;
end
chunks = chunks(~cellfun(@isempty, chunks));
if isempty(chunks)
    tblEvent = emptyTblEvent();
    return;
end
tblEvent = vertcat(chunks{:});

% Move sbjID to the front for readability.
varOrder = ['sbjID', setdiff(tblEvent.Properties.VariableNames, ...
    {'sbjID'}, 'stable')];
tblEvent = tblEvent(:, varOrder);
end


function tbl = emptyTblEvent()
tbl = table( ...
    categorical(strings(0, 1)), ...
    categorical(strings(0, 1), {'Cyto', 'Mito'}), ...
    zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    'VariableNames', ...
    {'sbjID', 'compartment', 'start', 'stop', 'amp', 'dur', 'int'});
end
