classdef reduct_displayer_IED_disp < multi_window_lock
% REDUCT_DISPLAYER_IED_DISP - A class to display the original IED data, in interaction with reduct_displayer.
% This class is a subclass of multi_window_lock.
%
% INPUTS: (Name-Value)
%       reduct_app: reduct_displayer - the reduct_displayer object that this class is connected to.
%       ied: IED.data - the IED data to display.
%       pos2start: double - the position in the recording to start from, in samples.
%       temp_events: cell - event types to display in addition to the detected events, that should not be saved in the ied object. See also IED_data2disp.
%       window_sizes: double - the size of the time windows for each data source, in seconds. The last element is the size of the zoom out window.
%       only_init_accepted: logical - if true, only events that are marked as "accepted" will be displayed. NOTE: UNCLEAR BEHAVIOR.

    properties
        reduct_app
    end
    properties
        last_pos
        current_pos
    end

    methods
        function app = reduct_displayer_IED_disp(options)
            arguments
                options.reduct_app reduct_displayer
                options.ied IED.data
                options.pos2start double {mustBeInteger, mustBePositive}
                
                options.temp_events
                options.window_sizes (:,1) double
                
                options.only_init_accepted (1,1) logical = false
            end

            % set deafult windows
            nData = height(options.ied.data_sources);
            if ~isfield(options, "window_sized")
                % default for everything is 1 sec, exept zoom out which is 10
                options.window_sizes = [ones(nData+1,1) ; 10];
            end

            % validate window_sizes
            nWins = numel(options.window_sizes);
            if nWins ~= (nData+2)
                error('Number of windows must match number of data+2 (for displaying zoom in & events)')
            end
            
            % transform window_sizes to init_time
            ied = options.ied;
            wins = options.window_sizes;
            current_pos = options.pos2start;
            detect_time = ied.pos(current_pos)./ied.fs;
            options.init_time = detect_time + [-wins, wins]./2;
            options = rmfield(options, ["window_sizes", "pos2start"]);
            
            % set aside the reduct_app
            reduct_app = options.reduct_app;
            options = rmfield(options, "reduct_app");
            % pass to super
            name_val = namedargs2cell(options);
            app@multi_window_lock(name_val{:})
            
            app.reduct_app = reduct_app;
            app.current_pos = current_pos;

        end

        function base_app_build(app, varargin)
            % build the parts of the app that are not dependent on input data
            
            %%%%%% Init app parts

            % deal with model
            if nargin == 2 && isa(varargin{1},"data2disp_model")
                app.model = varargin{1};
            else
                % solve weird bug - class do not init new data2disp_general_model upon start
                app.model = IED_data2disp();
            end
            app.data_loaded_listener = listener(app.model, "data_loaded", @(~,~) app.rebuild_app);

            % deal with views
            app.views.data_ax_view  = IED_data_ax(app, app.model);
            app.views.event_map     = event_map_time_win_fixed(app, app.model);

            app.controlers.time_window_controlers   = time_window_controlers_fixed(app, app.model);
            app.controlers.Y_lim_controlers         = Y_lim_controlers(app);
            app.controlers.events_controlers        = events_controlers(app, app.model, app.views.event_map);
            app.controlers.reduct_int_controlers    = reduct_int_controlers(app, app.model, app.views.event_map);

            %%%%%% Create Containers
            % create figure
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [1 1 1];
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'Event Marker';
            app.UIFigure.KeyReleaseFcn = @(~,evt) app.key_release_fcn(evt);
            app.create_containers()
            
            %%%%%% Set View Objects
            app.views.data_ax_view.set_view();
            app.views.event_map.set_view(...
                "local_map_loc", {app.containers.data_pannel, 1, 2}, ...
                "global_map_loc", {app.containers.controls_pannel, 1, [1 3]});
            
            %%%%%% Set controlers that do not depend on view changes
            app.controlers.time_window_controlers.set_controlers(...
                "time_jumper", {app.containers.controls_pannel 2 1}, ...
                "play_speed",  {app.containers.controls_pannel 2 2}, ...
                "play_button", {app.containers.controls_pannel 2 3},  ...
                "ROI_time", app.views.event_map.ROI_time);
            app.controlers.Y_lim_controlers.set_controlers(...
                "zoom_together_tree", {app.containers.data_pannel 1 1})
            app.controlers.events_controlers.set_controlers(...
                "mark_event_button",{app.containers.controls_pannel 2 5}, ...
                "event_chooser_dropdown",{app.containers.controls_pannel 1 5},...
                "marking_area",{app.containers.data_pannel, 1, 2}, ...
                "listen2line_click", true)
            app.controlers.reduct_int_controlers.set_controlers(...
                "return2curr_button", {app.containers.reduct_app_pannel 1 1}, ...
                "return2prev_button", {app.containers.reduct_app_pannel 2 1})
            

        end

        function create_containers(app)

            % add all super containers
            create_containers@multi_window_lock(app)
            
            % expend controls_pannel
            app.containers.controls_pannel.ColumnWidth = {'0.4x', '0.4x', '0.4x', ... % Area for time2knob
                '0.1x', '0.1x'}; % other buttons

            % add a sub-container for IED buttons
            reduct_app_pannel = uigridlayout(app.containers.controls_pannel);
            reduct_app_pannel.ColumnWidth = {'1x'};
            reduct_app_pannel.RowHeight = compose('%.4fx',[0.25,0.25]);
            reduct_app_pannel.Layout.Row = [1 2];
            reduct_app_pannel.Layout.Column = 4;
            reduct_app_pannel.BackgroundColor = [1 1 1];
            reduct_app_pannel.Padding = [0 0 0 0];
            app.containers.reduct_app_pannel = reduct_app_pannel;
        end

        function rebuild_app(app)
            % add protection to temporary types (i.e. clusters & labels)
            rebuild_app@multi_window_lock(app)

            app.controlers.events_controlers.set_controlers(...
                "protected_event_types", ["Detected";app.model.temp_event_types])
        end
    end
    
    %%%%%%%%%%%%% Helper Functions
    methods (Hidden)
        function [data_ax_p, row_loc, col_loc, YLim_editboxes_col_loc] = create_new_data_axes_pos(app,nLabels)
            % create where each data axe should go. Assume that data_axes
            % are in the same container as event_ax, and they share a
            % column. Make sure each axe have their own row.
            % Make sure that "Zoom Out" is at the bottom row compared to all other signals.

            % collect everything the same as parent
            [data_ax_p, row_loc, col_loc, YLim_editboxes_col_loc] = create_new_data_axes_pos@multi_window_lock(app, nLabels);
            
            % find the zoom-out & zoom-out source
            data_labels = app.model.data_labels;
            zoom_out_loc = contains(data_labels, " Zoomed Out");
            zoomed_out_source = extractBefore(data_labels(zoom_out_loc)," Zoomed Out");
            zoomed_out_source_loc = contains(data_labels, zoomed_out_source) & ~zoom_out_loc;

            % reorder so the "zoom out" get the largest number, and its
            % source the secound largest
            sorted_rows = sort(vertcat(row_loc{:}));
            new_rows = nan(size(sorted_rows));
            new_rows(zoom_out_loc) = sorted_rows(end);
            new_rows(zoomed_out_source_loc) = sorted_rows(end-1);
            new_rows(~(zoomed_out_source_loc | zoom_out_loc)) = sorted_rows(1:(end-2));

            % return to cell
            row_loc = num2cell(new_rows);

        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Keyboard Responce %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function key_release_fcn(app, evt)
            % responce to key release
            key_pressed = string(join([evt.Modifier evt.Key],'+'));
            app.controlers.time_window_controlers.key_responce(key_pressed)
            app.controlers.events_controlers.key_responce(key_pressed)
            app.controlers.Y_lim_controlers.key_responce(key_pressed)
        end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% Set Get Methods %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        function set.current_pos(app, val)
            app.last_pos = app.current_pos;
            app.current_pos = val;
            app.model.move_time(app.model.ied.pos(app.current_pos)/app.model.ied.fs)

        end
    end

    
end