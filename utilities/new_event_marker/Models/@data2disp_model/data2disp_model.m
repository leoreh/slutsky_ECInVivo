classdef data2disp_model < handle
    
    properties (Abstract)
        % hold the data to display.
        % Is a table, with variables "data" (cell) & "fs" (double), and
        % row names are data labels. 
        % Therefore, indexable by data_sources.data{data_label},
        % and you can find all sources labels by data_sources.Properties.RowNames
        data_sources table
    end

    properties
        % time of each event.
        % dictionary, key is the label of the event, field data is the
        % time stamps (numeric, 2 X nEvents).
        % Time stamps of each event is 1 row.
        % Events of the same type can't overlap.
        events_time (1,1) cell_dict  = cell_dict("Event",{double.empty(0,2)})
    end

    properties (Dependent)
        longest_data % find time in sec of the longest of the data_sources
        data_labels  % collect all data labels as string row vec
        events_types % collect all event types as string row vec
    end
    
    properties
        % current data to display
        % is a table, with variables "data", "actual_time_window", "req_time_window",
        % and rows are data labels. Therefore, indexable by data2disp.data{data_label}.
        %   data - actual channel data, numeric vector.
        %   actual_time_window - time window that matches the data, 2
        %                        element vector. Should be always wider
        %                        than "req_time_window".
        %   req_time_window - 
        data2disp table = table({},{},{},'VariableNames',["data", "actual_time_window", "req_time_window"])
        
        % same as events_time, but contains only the events in events_time2disp
        events2disp   (1,1) cell_dict = cell_dict("Event",{double.empty(0,2)}) 
        events_time2disp (1,2) double = [nan nan] % start & end time of current window for events view
    end

    events
        in_win_change % some data inside the view window may have changed
        data_loaded   % the source data in the object may have changed
        events_changed % when events data were changed (event added / removed)
        events_types_changed % when event types were changed (a whole type was added / removed)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% Deal With Loading New Data %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % activating "data_loaded" event
    methods
        function add_data(obj, data, options)
            % Adds data to the data2disp_model object and jumps to the first second.
            %
            % Syntax:
            %   add_data(obj, data, options)
            %
            % Input Arguments: (Positional)
            %   - obj: The data2disp_model object.
            %   - data: (repeating) The data to be added. Can be numeric matrix, file path, or memmap_row.
            %           If a single data is not a vector, it'll be separated by rows.
            %           If file path is given, file_map must be given as well, in a matching cell.
            %
            % Optional Arguments: (Name-Value Pair)
            %   - labels: Cell array, contining cellstr vectors. Each element corresponds to a data source. Inside each cell, element corresponds to a data row.
            %             If a cell is empty, default labels are used. Default repmat({''},1,numel(data)).
            %         For example, `{'vector_label', {'matrix_row1_label', 'matrix_row2_label', 'matrix_row3_label'}}`.
            %   - fs: Sampling frequency of each data. Default is 1250.
            %   - file_map: Cell array of file paths for each data. Default is an empty cell array.
            %   - events_time: Struct or dictionary containing event information. Default is an empty struct.
            %   - rmv_sets: Cell array of sets to be removed. Default is an empty cell array.
            %   - init_time: numeric (nData+1)-by-2, [start stop] time for
            %                model initation for each data source and events.
            %                Can also be a 1-by-2 vector, so same window
            %                will be used by all sources.
            %                
            %
            % Notes:
            %   - If no data is given, a dialog is opened to collect the data.
            %   - If data is given, it replaces any existing data in the object.
            %   - The function also sets the data2display to the first second and notifies listeners about the data loading and events changes.

            %%%%%%%%%%% Input declaration %%%%%%%%%%%%
            arguments
                obj
            end
            arguments (Repeating)
                data
            end
            arguments
                options.labels      (1,:) cell   = repmat({''},1,numel(data))
                options.fs (1,:) double {mustBeNumeric,mustBePositive} = ones(1,numel(data))*1250
                options.file_map    (1,:) cell   = cell(1,numel(data))
                options.events_time (1,:) {mustBeA(options.events_time,["struct","cell_dict", "dictionary"])} = struct("Event",[])
                options.rmv_sets    (1,:) cell   = cell(1,numel(data))
                options.init_time = [0 1];
            end

            % collect & validate source data
            if isempty(data) || ((numel(data) == 1) && isempty(data{:}))
                % no data given, open dialog. Note that this will only add
                % data (or events) to object, if data exists, and do not
                % replace.

                % collect data
                [data, options_from_dialog] = obj.dialog_add_data();
                
                % do not make any changes if user aborted
                if isempty(data)
                    return
                end

                % save data, everything is already validated
                obj.data_sources = data;
                obj.events_time = options_from_dialog.events_time;
            else
                % extract data from what was given. Note that this can replace existing data
                obj.data_sources = obj.extract_validate_data(data, options);
                obj.events_time = obj.validate_events(options.events_time);
            end
            
            % preset data2display to the first second, 
            % do not notify about it to maintain notify order
            obj.move_time(options.init_time, [], false)

            % finish loading, notify about it by order:
            notify(obj, "data_loaded")
            notify(obj, "events_types_changed")
            % notify(obj, "events_changed")
            notify(obj, "in_win_change") % note that this include any "events_changed" responce
        end
        
    end

    % helper functions (do not notify)
    methods (Hidden)
        [data, events_struct, options] = dialog_add_data(obj);

        function [extracted_data] = extract_validate_data(~, data, options)
            % validate that input data is OK.
            % INPUTS:
            %   data    - cell array, each cell is the input data - numeric
            %             mat, file path or memmap_row.
            %             If a single data is not a vector, it'll be
            %             separated by rows.
            %   options - struct, with fields:
            %       labels   - cell array of cellstr that match data, labels for each
            %                  input data. Each subcell match a data
            %                  channel (row).
            %       fs       - numeric array that match data, sampling freq
            %                  for each input data.
            %       rmv_sets - cell array that match data, in each cell
            %                  what rows to ignore - leave empty if to
            %                  ignore non.
            %       file_map - cell array that match data, in each cell how
            %                  the corresponding file is mapped (cell array).
            %                  See memmapfile 'Format' option with 3 cell,
            %                  or map_creator.
            %                  Required when matching "data" is file.
            % OUTPUT:
            %   extracted_data - 
            %             table containing data sources, rows names are
            %             labels, variables are "data" (cell) & "fs" (double).

            %%%%%%%%%%%% Input validation %%%%%%%%%%%%
            % make sure that options match data
            nData = numel(data);
            if nData == 0
                % just return empty upon no data entry
                extracted_data = table('Size',[0 2],'VariableTypes',["cell","double"],'VariableNames',["data","fs"]);
                return
            elseif numel(options.labels) ~= nData
                error('mismatch number of labels & data')
            elseif numel(options.fs) ~= nData
                error('mismatch number of fs & data')
            elseif numel(options.file_map) ~= nData
                error('mismatch number of file_formats & data')
            elseif numel(options.rmv_sets) ~= nData
                error('mismatch number of rmv_sets & data')
            end

            %%%%%%%%%%%% Extract data & data labels %%%%%%%%%%%%
            for iData = nData:-1:1
                current_data = data{iData};


                %%%% extract data

                % treat numeric & file separately
                if isa(current_data,'memmap_row')

                    % memmap_row is already ordered as needed -
                    % just save it in a cell
                    data{iData} = {current_data};

                elseif isnumeric(current_data)

                    % error if you can't seperate
                    if ~ismatrix(current_data)
                        error('Data input must be vector or 2d matrix, so each row will be different data set. Input %d does not match that',iData)
                    end

                    % seperate rows
                    if isvector(current_data) && ~isrow(current_data)
                        % force row vector
                        current_data = current_data(:)';
                    end
                    data{iData} = num2cell(current_data,2)';

                elseif (ischar(current_data) || isstring(current_data) || iscellstr(current_data)) && isfile(current_data)

                    % convert to maped file
                    maped_file = memmapfile(current_data,"Format",[options.file_map{iData} 'Mapped']);

                    % error if you can't seperate
                    if ~ismatrix(maped_file.Data.Mapped)
                        error('Data input must be vector or 2d matrix, so each row will be different data set. Input %d does not match that',iData)
                    end

                    % seperate rows
                    for iRow = size(maped_file.Data.Mapped,1):-1:1
                        %                         row_ref{iRow} = @(indx) maped_file.Data.Mapped(iRow,indx);
                        row_ref{iRow} = memmap_row(maped_file,iRow);
                    end
                    data{iData} = row_ref;

                else
                    error("Data %d is of unsupported type (should be numeric or binary file)",iData)
                end

                % remove unwanted datasets -
                % delete cells that match the rows that were separated earlier.
                data{iData}(options.rmv_sets{iData}) = [];


                %%%% extract labels
                % create labels if needed
                if isempty(options.labels{iData})
                    options.labels{iData} = {sprintf('data%d',iData)};
                end

                % if label is not cellstr in cell but something else that
                % can be used - user want this label to be the base to all
                % the sub-channels, or we have only 1 channel. 
                % Wrap it as cellstr in cell so expending will work as planned
                % later.
                if ~iscellstr(options.labels{iData})
                    options.labels{iData} = cellstr(string(options.labels{iData}));
                end

                % expend to match separation if needed
                if numel(data{iData}) ~= numel(options.labels{iData})
                    if numel(options.labels{iData}) == 1
                        % 1 label for all - expend it by adding row number
                        options.labels{iData} = cellstr((options.labels{iData} + "_" + string(1:numel(data{iData})) ));
                    else
                        error('Mismatch nLabels & data after row seperation, in data %d',iData)
                    end
                end
            end
            
            % make sure label is not "Events" - reserved for events
            if ismember("Events",options.labels{iData})
                error('Label "Events" is reserved for events, and can not be used as data label, in data %d',iData)
            end

            % place data and labels in app, so each cell is a different data source
            data_sources = [data{:}];
            data_labels  = [options.labels{:}];
            data_fs      = repelem(options.fs,cellfun(@numel,data));

            %%%% validate data length
            % get length of each data, validate they match (within 0.5 sec error)
            data_time = cellfun(@(x,fs) length(x)./fs,data_sources,num2cell(data_fs)); % time in sec
            longest_data = max(data_time);
            if any(~ismembertol(data_time,longest_data,0.5,"DataScale",1))
                % should be within 0.5 sec from ending
                error('Data time should match within 0.5 sec from longest')
            end

            % create data table
            extracted_data = table(data_sources', data_fs', 'VariableNames',["data","fs"],'RowNames',data_labels);
        end

        function [events_time] = validate_events(obj, events_time, merge_flag)
            % validate & collect the events.
            %
            % INPUTS:
            %   app           - event_marker app.
            %   events_time   - struct, each field is an event type OR
            %                   dictionary / cell_dict, each key is event type.
            %                   each type contain a matrix with 2 cols,
            %                   col 1 is event start col 2 is event end.
            %   merge_flag    - logical flag, if true will merge
            %                   overlapping events to create single long
            %                   events.
            %                   Defualt true;
            %
            % OUTPUTS:
            %   events_time   - simillar to input, but always cell_dist.
            %
            % The function will error if validatation is failed.
            
            if ~exist("merge_flag","var")
                merge_flag = false;
            end
            % convert input to cell_dict, if needed
            if isa(events_time, "struct")
                new_keys = fieldnames(events_time);
                new_values = struct2cell(events_time);
                events_time = cell_dict(string(new_keys), new_values);
            elseif isa(events_time, "dictionary") || isa(events_time, "containers.Map")
                events_time = cell_dict(events_time);
            end
            % collect info
            event_types = keys(events_time);
            event_types = string(event_types(:)');
            for iType = event_types
                if isempty(events_time{iType})
                    % make sure it is an empty vector with 2 columns
                    events_time{iType} = double.empty(0,2);
                    continue
                end

                % validate no event start before it ends
                if any(events_time{iType}(:,1) > events_time{iType}(:,2))
                    error("Event type %s have an event start (col 1) is after event stop (col 2)",iType)
                end
                
                if merge_flag
                    % merge overlapping events
                    events_time{iType} = merge_events(events_time{iType});
                end

                % sort events by start time
                events_time{iType} = sortrows(events_time{iType},[1 2],"ascend");

                % validate no event is placed after data ends ./ before it starts
                if any(obj.longest_data < events_time{iType}(end,2)) || any(events_time{iType}(1,1) < 0)
                    error('Event type %s have at least 1 event is after end of data or before 0', iType)
                end
            end
        end
    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% Deal With Moving to New Time %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % activating "in_win_change" event
    methods
        function move_time(obj, new_time, Sources, notify_flag)
            % collect data in new window, while maintaining the same window center. for all sources.
            % Note: window will be trimmed when it is out of bond, 
            %       see time2ind.
            % INPUTS:
            %   obj - caller object
            %   new_time - 2 options:
            %              a) numeric scalar, window center to move all sources to.
            %              b) numeric nSources X 2, window start-stop times (Sorted!).
            %                 data_sources in Sources will be moved to the matching window, 
            %                 all others will only change center.
            %                 All windows must have the same center (mean).
            %                 If new time is a 1 by 2, it will be used for
            %                 all sources given.
            %   Sources - cellstr or strings, data_sources labels to update time window, see new_time.
            %             Could also include "Events" to change events2disp, otherwise it will also change center.
            %             default is all "Events" & than all data_sources (in that order).
            %             Note that this has no effect if new_time is a scalar.
            %   notify_flag - logical flag, if to notify in_win_change event.
            %                 default true.
            
            %%%%%%%% Manage Input
            % in no source was specified, return all of them
            if ~exist("Sources","var") || isempty(Sources)
                Sources = ["Events" obj.data_labels];
            end
            
            % default notify_flag to true
            if ~exist("notify_flag","var")
                notify_flag = true;
            end

            % make sure Sources are row vector
            Sources = Sources(:)';

            % error if any source do not exist
            sources_not_exists = ~ismember(Sources, [obj.data_labels "Events"]);
            if any(sources_not_exists)
                 error('data2disp_model:move_time:BadSrcReq',...
                     'Sources {%s} does not exist.\nTo add data, use add_data function', ...
                     string(join(Sources(sources_not_exists),', ')))
            end
            
            % validate new_time
            new_time_center = mean(new_time,2);
            if ~issorted(new_time,2)
                error('data2disp_model:move_time:BadTimeReq','new_time must be sorted')
            elseif any(~ismembertol(new_time_center, new_time_center(1)))
                error('data2disp_model:move_time:BadTimeReq','All windows must have the same center')
            end
            new_time_center = new_time_center(1);

            % match number of new_time & sources
            nSources = numel(Sources);
            if all(size(new_time) == [1 2]) || isscalar(new_time)
                new_time = repmat(new_time,[nSources, 1]);
            elseif size(new_time,1) ~= nSources
                error('data2disp_model:move_time:BadTimeReq',...
                    'new time must be a scalar, a 2 elements row vec, or have 1 row each requested source')
            end
            
            %%%%%%%% Move time windows
            % collect data matching the new time window we moved to
            warning('off','MATLAB:table:RowsAddedExistingVars') % do not warn on new row names - there is data in data_sources that isn't yet in data2disp
            for iSource = obj.data_labels
                % current_new_time = new_time{iSource};
                
                % seperate between sources that were asked to move to new time, and those that were asked to move to new center
                req_source_check = ismember(Sources, iSource);
                if any(req_source_check)
                    current_new_time = obj.center2window(new_time(req_source_check,:), iSource);
                else
                    current_new_time = obj.center2window(new_time_center, iSource);
                end

                % transfer from time to ind for that dataset
                [dataset_ind, actual_new_time] = obj.time2ind(iSource, current_new_time);
                
                % register new window
                obj.data2disp.data{iSource} = obj.data_sources.data{iSource}(dataset_ind);
                obj.data2disp.actual_time_window{iSource} = actual_new_time;
                obj.data2disp.req_time_window{iSource} = current_new_time;
            end
            warning('on','MATLAB:table:RowsAddedExistingVars')
            
            % update events window
            req_source_check = ismember(Sources, "Events");
            if any(req_source_check)
                current_new_time = obj.center2window(new_time(req_source_check,:), "Events");
            else
                current_new_time = obj.center2window(new_time_center, iSource);
            end
            % deal with scalar new time (window center), save new time
            obj.events_time2disp = obj.center2window(current_new_time, "Events");

            % run through all types, find matching events
            for iType = obj.events_types
                obj.events2disp{iType} = obj.find_event2disp(iType);
            end

            % finished changing, notify about it, unless user asked not to
            % (initiating model)
            if notify_flag
                notify(obj,'in_win_change')
            end
        end    
    
        function add_event(obj, event_t, event_type)
            % take new event time, and write it down to the right event
            % type.
            % In case of overlaps, new event replace old ones.
            %
            % event_t    - numeric 2 elem, [start stop].
            % event_type - string scalar

            % if event type do not exist, add it
            if ~ismember(event_type, obj.events_types)
                obj.add_event_type(event_type)
            end

            % collect events that match this event
            relv_events = obj.events_time{event_type};

            % find the relevant event
            % check if new event time overlap with existing event % saved for legacy
            %{
            events_starts = relv_events(:,1);
            events_ends   = relv_events(:,2);

            % overlapping cases: 
            % new event contain other event, fully or pariatly:
            %   any of old event time, is between new event time
            contain_overlap_start   = (event_t(1) <= events_starts) & (events_starts <= event_t(2));
            contain_overlap_end     = (event_t(1) <= events_ends)   & (events_ends <= event_t(2));
            new_event_contain_old = contain_overlap_start | contain_overlap_end;
            % new event is contained fully in old event:
            %   start of new event time is between old event times
            %   (could be end of new event, same thing for fully contained case)
            old_event_contain_new = (events_starts <= event_t(1)) & (event_t(1) <= events_ends);

            events2replace = new_event_contain_old | old_event_contain_new;
            %}
            
            % find overlaps and replace them -
            % a better solution could be to only calculate for 1 event, but
            % the general solution is fast enough & vectorized
            overlap_mat = check_overlap([relv_events; event_t]);
            events2replace = overlap_mat(:,end);
            if any(events2replace)
                % if overlap, replace the existing events instead of insert anew
                relv_events(events2replace,:) = [];
            end

            % insert new event into event list
            obj.events_time{event_type} = sortrows([relv_events; event_t],[1 2],'ascend');

            % see if event should be added to current window
            obj.events2disp{event_type} = obj.find_event2disp(event_type);

            % notify that event list was changed
            notify(obj, 'events_changed')
        end
    
        function rmv_event(obj, time_point, event_type)
            % remove all events of a type that contain a specific
            % time point.
            % contained_t - numeric scalar, time point to check for intersections.
            % event_type  - string scalar

            [~, relv_event_log] = point2event(obj, time_point, event_type);

            % insert events without deleted
            obj.events_time{event_type}(relv_event_log,:) = [];

            % Check if this was tha last event 
            if isempty(obj.events_time{event_type})
                % maintain 2 cols interface
                obj.events_time{event_type} = double.empty(0,2);
            else
                % make sure events are ordered
                obj.events_time{event_type} = sortrows(obj.events_time{event_type},[1 2],"ascend");
            end

            % see if event should be removed in current window
            obj.events2disp{event_type} = obj.find_event2disp(event_type);

            % notify that event list was changed
            notify(obj, 'events_changed')
        end
        
        function add_event_type(obj, event_type)
            % add a new event
            %   event_type - string scalar

            if ~ismember(event_type, obj.events_types)
                obj.events_time{event_type}     = double.empty(0,2);
                obj.events2disp{event_type}  = double.empty(0,2);
            else
                error("data2disp_model:add_event_type:event_exists",...
                    "Event %s already exist in the data", event_type)
            end

            % notify that event list was changed
            notify(obj, 'events_types_changed')
        end

        function rmv_event_type(obj, event_type)
            % remove all events of a certain type.
            %   event_type - string scalar

            if ismember(event_type, obj.events_types)
                obj.events_time(event_type)     = [];
                obj.events2disp(event_type)  = [];
            else
                error("data2disp_model:rmv_event_type:no_event_exists",...
                    "Event %s do not exist in the data", event_type)
            end

            % notify that event list was changed
            notify(obj, 'events_types_changed')
        end
    
    end
    
    
    % non-hidden helper methods (do not notify)
    methods
        function [new_time] = center2window(obj, win_center, iData)
            % if given a scalar win_center, expend so window size will be
            % kept but centered on a new point.
            % if win_center isn't scalar, return it as it is.
            % note: using center will fail during initation. if data2disp.req_time_window &
            %   events_time2disp are empty, you need to use a full time
            %   window.

            if isscalar(win_center)

                % collect current window
                if strcmpi(iData,"Events")
                    view_window = obj.events_time2disp;
                else
                    view_window = obj.data2disp.req_time_window{iData};
                end

                % create start-stop times with same center & window size
                win_size = diff(view_window);
                new_time = win_center + [-win_size win_size]./2;
            else
                % time is already start-stop, force it to row vec
                new_time = win_center(:)';
            end
        end
        
        function [relv_event, relv_event_log] = point2event(obj, time_point, event_type)
            % check if a single time point overlaps with any event of
            % certain type, if so retrun it. Else return empty.

            events_starts = obj.events_time{event_type}(:,1);
            events_ends   = obj.events_time{event_type}(:,2);
            relv_event_log = (events_starts <= time_point)  & (time_point <= events_ends);
            relv_event = [events_starts(relv_event_log),events_ends(relv_event_log)];
        end
    
    end
    
    % helper functions (do not notify)
    methods (Hidden)
        function [dataset_ind, actual_new_time] = time2ind(obj, iSource, new_time)
            % Find the samples corresponding to the specified time window (new_time) for the specified dataset (iSource).
            % Adjust the window to include whole sample indices if necessary.
            % Add samples to match the window size if the number of samples is insufficient.
            % Trim the signal if the requested samples are out of range.
            % Note: new_time must be in ascending order.
            
            % collect sample window
            fs = obj.data_sources.fs(iSource);
            ind_base = new_time .* fs;
            ind_base = [floor(ind_base(1)), ceil(ind_base(2))];
            
            % find the mismatch between the requested time and the actual time
            actual_new_time = ind_base/fs;
            if any(actual_new_time ~= new_time)
                % find the difference
                time_diff = actual_new_time - new_time;
                % find the number of samples to add
                sample_diff = time_diff * fs;

                % sample diff is a 2x1 array, with how many samples to add to each side.
                % add the samples.
                ind_base(1) = ind_base(1) - ceil(sample_diff(1));
                ind_base(2) = ind_base(2) + ceil(sample_diff(2));
            end

            % deal with edge cases - shrink window as needed.
            data_len = length(obj.data_sources.data{iSource});
            % end after dataset size
            if ind_base(2) > data_len
                % trim right
                ind_base(2) = data_len;
            end
            % starting before time 0
            if ind_base(1) < 1
                % trim left
                ind_base(1) = 1;
            end
            
            % recalculate actual_new_time after all the changes
            actual_new_time = ind_base/fs;

            % create ind
            dataset_ind = ind_base(1):ind_base(2);
        end
        
        function events2disp = find_event2disp(obj, event_type)
            % check which events of certain event type are currently in
            % display window
            
            % collect all events of type
            all_event_t = obj.events_time{event_type};

            % collect what events occure in the time frame
            events2exclude = [all_event_t(:,2) < obj.events_time2disp(1), all_event_t(:,1) > obj.events_time2disp(2)];
            events2exclude = any(events2exclude,2);
            events2disp = all_event_t(~events2exclude,:);
        end
    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% Get Set Methods %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function res = get.longest_data(obj)
            % calculate the length in sec of the longest data source

            data_time = cellfun(@(x,t) length(x)./t, obj.data_sources.data, num2cell(obj.data_sources.fs)); % time in sec
            res = max(data_time);
        end
    
        function res = get.data_labels(obj)
            % collect data labels, force them to be string row vec
            
            res = obj.data_sources.Properties.RowNames;
            res = string(res(:)');
        end
    
        function res = get.events_types(obj)
            % collect all events types, return them as a string row vec

            res = keys(obj.events_time);
            res = string(res(:)');
        end
    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% Save & Load %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function save_only_evt_data(obj, save_loc)
            % save the event's data extracted from the model

            % source_data = obj.data_sources.data;
            fs = obj.data_sources.fs;
            datas_label = obj.data_labels;
            event_times = obj.events_time;
            events_data = obj.collect_evt_data();

            save(save_loc,...
                'event_times','events_data','datas_label','fs')
        end
    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% Utilites %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        function [new_time] = refit_window(obj, new_time)
            % make sure time is limited by data limits.
            %   INPUTS:
            %       obj - caller model object.
            %       new_time - [start stop] time window.
            %
            %   OUTPUT:
            %       new_time - same time window, moved so window does not
            %                  pass longest data or 0. Keeping window size
            %                  constant if possible.
            %
            % NOTE1: This function provides convenience for model users by
            %        adjusting the time window, which differs from the basic
            %        model behavior, that trims the window as needed.
            % NOTE2: This function isn't source specifc, but compares
            %        time to longest data & events_time2disp. 
            %        Inspect data2disp.time_window to consider specific
            %        data_source limits.

            t_end = obj.longest_data;
            win_size = diff(new_time);
            if diff(new_time) > t_end
                new_time = [0 t_end];
            end
            if new_time(2) > t_end
                new_time = [t_end-win_size, t_end];
            end
            if new_time(1) < 0
                new_time = [0, win_size];
            end
        end
    
        function [collected_data] = collect_data(obj, iSource, time2collect)
            % collect some time point without changing the window
            
            ind = obj.time2ind(iSource, time2collect);
            collected_data = obj.data_sources.data{iSource}(ind);
        end
        
        function [evt_data] = collect_evt_data(obj, sources2extract)
            % collect the data in all datasources tha matches events time
            % OUTPUT:
            %   events_data - struct:
            %       fields are event types, data is cell
            %       matrix - each row is different data_source, each col is a
            %       different event. Sources skipped will result in a cell row
            %       of emptys.
            
            % default sources to all data_labels
            if ~exist("Sources","var") || isempty(sources2extract)
                sources2extract = obj.data_labels;
            end
            
            % error if any source do not exist
            sources_not_exists = ~ismember(sources2extract, obj.data_labels);
            if any(sources_not_exists)
                 error('data2disp_model:collect_evt_data:BadSrcReq',...
                     'Sources {%s} does not exist.\nTo add data, use add_data function', ...
                     string(join(Sources(sources_not_exists),', ')))
            end

            % loop through event types, event & sources to fill evt_data
            % Note that MATLAB is column order, this is slightly faster
            % than traversing along the rows

            for iType = obj.events_types
                nEvents = size(obj.events_time{iType},1);
                
                for iEvent = nEvents:-1:1
                    event_t = obj.events_time{iType}(iEvent,:);
                    
                    for iSource = numel(obj.data_labels):-1:1
                        current_source = obj.data_labels(iSource);

                        % only fill sources that were requested - all else
                        % will remain empty rows
                        if ismember(current_source, sources2extract)
                            evt_data.(iType){iSource,iEvent} = obj.collect_data(current_source, event_t);
                        end
                    end

                end

            end
        end

    end

    methods (Static)
        function [event_vec] = events_nan_vec(events_time)
            % convert a list of events time, which is orginazed as nEvents
            % by 2 [start stop] times, into 1 vector. Each event will be
            % separated by nan from the next event.
            %
            % INPUT:
            %   events_time - nEvents X [start stop] numeric, each event
            %                 time.
            %
            % OUTPUT:
            %   event_vec - (3*nEvents) X 1 vector, all the events
            %               separated by nans.
            
            % add nan col
            nEvents = size(events_time,1);
            events_time = [events_time nan(nEvents,1)];

            % transpose & convert to vector.
            % transpose as MATLAB collection is column oriented.
            events_time = events_time';
            event_vec = events_time(:);
        end
    
    end
end