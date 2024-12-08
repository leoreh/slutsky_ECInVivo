function [data, options] = dialog_add_data(obj)
    % interactivly add data to model using data_loader app

    loader_app = data_loader(...
        'non_editable', 'loadable', @(varargin) loader_validate_input(obj, varargin{:}), ...
        'editable', 'Role', categorical(string.empty,["Signal","Events"],'Protected',true), ...
        'editable', 'Data_Label', string.empty, ... % use spaces 2 seperate cells
        'editable', 'fs', [], ...
        'editable', 'rmv_sets', string.empty, ... % use spaces 2 seperate cells
        'editable', 'File_Map_precision', categorical(string.empty,["int8","int16","int32","int64","uint8","uint16","uint32","uint64","single","double"],'Protected',true), ...
        'editable', 'File_Map_nChannels', [], ...
        'non_editable', 'signal_length', @loader_sig_duration, ...
        'non_editable', 'event_minmax', @loader_event_minmax, ...
        'file_filters',{'*.dat' '*.lfp'}, 'workspace_filters', ...
        @(x) ismember(x.class,["struct","dictionary", "int8","int16","int32","int64","uint8","uint16","uint32","uint64","single","double","memmap_row"])...
        );
    % loader_app.UIFigure.Position(3) = 1519;
    movegui(loader_app.UIFigure,'center')

    % make sure to close loader_app if exit on any reason
    cl2 = onCleanup(@() delete(loader_app));

    % loop until table is validated
    while true

        % wait until window is ready to be closed - closing call uicontinue
        uiwait(loader_app.UIFigure)

        % collect inputs 2 load, validate it
        load_table = loader_app.UITable.Data;
        if isempty(load_table)
            % user did not ask to load anything, abort
            [data, options] = deal([]);
            return
        else
            % validate and parse inputs 2 load
            [res, data, options] = loader_validate_input(obj,load_table);
            if res == "✓"
                % we can load, return
                break
            else
                % raise alert, run loop again
                uialert(loader_app.UIFigure,["Can't load, due to:" res "There may be more issues, look at the whole table!"],"Load Failed","Icon","error","Modal",true)
            end
        end
    end
end


function res = loader_sig_duration(row_in)
    % function for info col "signal_length" at "data_loader":
    % calc signal length in time, return duration scalar

    % answer for no input
    if nargin == 0
        res = duration.empty;
        res.Format = "hh:mm:ss.SSS";
        return
    end

    % validate input can be calculated
    if row_in.Role ~= "Signal"
        error('duration is only relevant to signals')
    end
    if ~(row_in.fs > 0)
        error('fs must be non negative')
    end

    % read by type
    if row_in.Source == "Workspace"
        % input is workspace var

        % collect input size
        % Before calling the whos function from evalin, we'll do a
        % security check that the whos function is the one provided
        % by Matlab.  This reduces the risk of using evalin().
        if isempty(regexp(which('whos'), 'built-in *.', 'once'))
            error('The built-in whos() function is being overshaddowed.  Do not proceed.')
        end
        baseVarInfo = evalin('base', sprintf('whos("%s")',row_in.Name));

        % get samples per channel
        if numel(baseVarInfo.size) ~= 2
            error('Input var must be vector or matrix')
        end
        if any(baseVarInfo.size == 1)
            % var is a vector, take length
            nSamps = max(baseVarInfo.size);
        else
            % var is matrix, will be seperated by rows
            nSamps = baseVarInfo.size(2);
        end

    else
        % input is file

        % validate it has mapping
        if ismissing(row_in.File_Map_precision) || ismissing(row_in.File_Map_nChannels)
            error('Missing file map to prase file')
        end

        % collect full file path
        file_path = fullfile(row_in.Source, join(row_in{:,["Name","Type"]},''));

        % collect samples per channel
        map_format = map_creator(file_path, row_in.File_Map_nChannels, row_in.File_Map_precision);
        nSamps = map_format{2}(2);
    end

    % calculate duration
    res = seconds(nSamps./row_in.fs);
    res.Format = "hh:mm:ss.SSS";

end

function res = loader_event_minmax(row_in)
    % function for info col "signal_length" at "data_loader":
    % find the time of the first & last event in the struct,
    % return duration 2 elements

    % answer for no input
    if nargin == 0
        res = duration.empty(0,2);
        return
    end

    % only work with events
    if row_in.Role ~= "Events"
        error('event_minmax only relevant to events')
    end

    % work only with workspace var & struct
    if row_in.Source ~= "Workspace" || row_in.Type ~= "struct"
        error('Only supported workspace loaded structs')
    end

    % collect the variable
    base_var = evalin('base',row_in.Name);

    % collect min & max for every event type
    events_min = structfun(@(x) min(x(:,1)),base_var);
    events_max = structfun(@(x) max(x(:,2)),base_var);

    % calculate duration
    res = seconds([min(events_min), max(events_max)]);
    res.Format = "hh:mm:ss.SSS";
end

function [res, data, options] = loader_validate_input(obj,table_in)
    % function for info col "signal_length" at "data_loader" & for
    % validation upon "data_loader" closing.
    % check if for a given data set in app & for specific input
    % table from "data_loader", we can load the data in the table.
    %
    % If we can load, return res "✓", else return the error text.
    %
    % If we can load, return also data, events & options that can
    %   be inserted into app.
    %
    % Use app.extract_validate_data & app.validate_events.

    % answer for no input
    if nargin == 1
        res = string.empty;
        [data, options] = deal([]);
        return
    end

    % seperate incoming table by role
    only_signals = table_in(table_in.Role == "Signal",:);
    only_events  = table_in(table_in.Role == "Events",:);
    only_undefined = table_in(ismissing(table_in.Role),:);

    % return if any undefined
    if ~isempty(only_undefined)
        bad_names = "{" + join(only_undefined.Name,', ') + "}";
        res = sprintf("%s: missing role",bad_names);
        [data, options] = deal([]);
        return
    end

    %%%%%% Validate Signals

    % collect signals input
    data = cell(1,height(only_signals));
    options.file_map  = cell(1, height(only_signals));
    options.rmv_sets   = cell(1, height(only_signals));
    options.labels     = cell(1, height(only_signals));
    for iSig = 1:height(only_signals)

        % parse the signal & accompaning data
        if only_signals.Source(iSig) == "Workspace"
            % collect from workspace.
            % note security issue if workspace is changed.
            data{iSig} = evalin("base",only_signals.Name(iSig));
        else
            % collect file info.
            data{iSig} = fullfile(only_signals.Source(iSig), join(only_signals{iSig,["Name","Type"]},''));

            % collect map for file
            if mod(only_signals.File_Map_nChannels(iSig),1) ~= 0 || ismissing(only_signals.File_Map_precision(iSig))
                bad_names = only_signals.Name(iSig);
                res = sprintf("%s: missing/ bad data for file map",bad_names);
                [data, options] = deal([]);
                return
            end
            options.file_map{iSig} = map_creator(data{iSig}, only_signals.File_Map_nChannels(iSig), only_signals.File_Map_precision(iSig));
        end

        % parse rmv_sets
        if ismissing(only_signals.rmv_sets(iSig))
            options.rmv_sets{iSig} = [];
        else
            options.rmv_sets{iSig} = split(strtrim(only_signals.rmv_sets{iSig}),[", ", " "])';
            options.rmv_sets{iSig} = str2double(options.rmv_sets{iSig});
            if any(mod(options.rmv_sets{iSig},1) ~= 0)
                bad_names = only_signals.Name(iSig);
                res = sprintf("%s: rmv_sets parse failed",bad_names);
                [data, options] = deal([]);
                return
            end
        end

        % parse labels
        if ismissing(only_signals.Data_Label(iSig)) || only_signals.Data_Label(iSig) == ""
            options.labels{iSig} = {};
        else
            options.labels{iSig} = split(strtrim(only_signals.Data_Label{iSig}),[", ", " "])';
        end
    end

    % validate collect fs
    miss_fs = ismissing(only_signals.fs) | (only_signals.fs <= 0);
    if any(miss_fs)
        bad_names = "{" + join(only_signals.Name(miss_fs),', ') + "}";
        res = sprintf("%s: missing fs",bad_names);
        [data, options] = deal([]);
        return
    end
    options.fs = only_signals.fs;

    % add app existing data
    data = [obj.data_sources.data ; data]';
    options.labels = [cellstr(obj.data_labels), options.labels];
    options.fs = [obj.data_sources.fs ; options.fs]';
    options.file_map = [cell(1,height(obj.data_sources)), options.file_map];
    options.rmv_sets = [cell(1,height(obj.data_sources)), options.rmv_sets];
    
    % validate & collect
    try
        [data] = obj.extract_validate_data(data, options);
    catch err
        res = string(err.message);
        [data, options] = deal([]);
        return
    end

    %%%%%% Validate events

    merged_events = obj.events_time;
    merge_entries = entries(merged_events);
    
    for iEvt_in = 1:height(only_events)
        % collect from workspace.
        % note security issue if workspace is changed.
        new_events = evalin("base",only_events.Name(iEvt_in));
        
        % combine events
        if isstruct(new_events)
            new_keys = fieldnames(new_events);
            new_values = struct2cell(new_events);
            merge_entries = [merge_entries; table(new_keys,new_values, 'VariableNames',["Key","Value"])];
        else
            merge_entries = [merge_entries; entries(new_events)];
        end
    
        % error on dup keys
        if height(merge_entries) ~= numel(unique(merge_entries.Key))
            res = ["Events in %s have duplicated events types with others." newline...
                "The duplicated types maybe already loaded, or in other events you want to load."];
            res = sprintf(res, only_events.Name(iEvt_in));
            [data, options] = deal([]);
            return
        end
    end
    % convert back to cell_dict
    merged_events = cell_dict(merge_entries.Key, merge_entries.Value);

    % validate events
    try
        options.events_time = obj.validate_events(merged_events);
    catch err
        res = string(err.message);
        [data, options] = deal([]);
        return
    end

    %%%%%% Mark everything is good
    res = "✓";

end

