classdef multi_window_lock < event_marker_gui
    % same as event_marker_gui, but have multiple diffrent time windows for
    % diffrent axes, all size-locked for the duration of the app.

    methods
        function app = multi_window_lock(varargin)
            % make sure user manually wrote what size is each window, if no
            % model was given
            
            if nargin == 1 && isa(varargin{1},"data2disp_model")
                % do nothing, superclass constructor will use it
            else
                % find init_time
                init_time_pos = cellfun(@(x) (isStringScalar(x) || (ischar(x) && isrow(x)))  && x == "init_time",...
                    varargin,'UniformOutput',true);
                if ~any(init_time_pos)
                    error('Must specify init_time')
                end
            end

            % call super constroctor
            app@event_marker_gui(varargin{:})
        end

        function base_app_build(app, varargin)
            % build the parts of the app that are not dependent on input data
            
            %%%%%% Init app parts

            % deal with model
            if nargin == 2 && isa(varargin{1},"data2disp_model")
                app.model = varargin{1};
            else
                % solve weird bug - class do not init new data2disp_general_model upon start
                app.model = data2disp_general_model();
            end
            app.data_loaded_listener = listener(app.model, "data_loaded", @(~,~) app.rebuild_app);

            % deal with views
            app.views.data_ax_view  = data_ax_view(app, app.model);
            app.views.event_map     = event_map_time_win_fixed(app, app.model);

            app.controlers.time_window_controlers   = time_window_controlers_fixed(app, app.model);
            app.controlers.Y_lim_controlers         = Y_lim_controlers(app);
            app.controlers.events_controlers        = events_controlers(app, app.model, app.views.event_map);
            app.controlers.general_controlers       = general_controlers(app, app.model);

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
                "global_map_loc", {app.containers.controls_pannel, 1, [1 4]});
            
            %%%%%% Set controlers that do not depend on view changes
            data_menus = uimenu(app.UIFigure,"Text","&Data");
            app.controlers.general_controlers.set_controlers(...
                "save_button", {app.containers.controls_pannel 2 4}, ...
                "add_close_req", true,...
                "add_data_menu", data_menus, ...
                "export_data_menu", data_menus, ...
                "save_loc_menu", data_menus)
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
                "marking_area",{app.containers.data_pannel, 1, 2})
        end

        function create_containers(app)
            % Create pannel_holder
            pannel_holder = uigridlayout(app.UIFigure);
            pannel_holder.ColumnWidth = {'1x'};
            pannel_holder.RowHeight = {'1x', '0.2x'};
            pannel_holder.BackgroundColor = [1 1 1];
            app.containers.pannel_holder = pannel_holder;

            % Create data_pannel
            data_pannel = uigridlayout(app.containers.pannel_holder);
            data_pannel.ColumnWidth = {'0.075x', '1x'};
            data_pannel.RowHeight = {'1x'};
            data_pannel.Layout.Row = 1;
            data_pannel.Layout.Column = 1;
            data_pannel.BackgroundColor = [1 1 1];
            app.containers.data_pannel = data_pannel;

            % Create controls_pannel
            controls_pannel = uigridlayout(app.containers.pannel_holder);
            controls_pannel.ColumnWidth = {'0.3x', '0.3x', '0.3x', '0.3x', ... % Area for time2knob
                '0.1x'}; % other buttons
            controls_pannel.RowHeight   = {'1x', '0.3x'};
            controls_pannel.Layout.Row = 2;
            controls_pannel.Layout.Column = 1;
            controls_pannel.BackgroundColor = [1 1 1];
            app.containers.controls_pannel = controls_pannel;
        end
    
        % data_loaded listener for: model
        function rebuild_app(app)

            %%%%%% remove existing axes & events
            data_ax_view = app.views.data_ax_view;
            data_ax_view.rmv_data_axes();
            app.views.event_map.rmv_map_lines();
            app.controlers.events_controlers.rmv_marking_area();
            % app.controlers.Y_lim_controlers.

            %%%%%% create uiprogressdlg if figure is visable
            if app.UIFigure.Visible
                progress_dlg = uiprogressdlg(app.UIFigure,"Indeterminate","on","Message","Loading data","Cancelable","off","Title","Loading data");
                cl = onCleanup(@() close(progress_dlg)); % remove progress_dlg block on any form of exit (i.e. error)
            end

            %%%%%% create new axes to match data
            data_labels = app.model.data_labels;
            nLabels = numel(data_labels);
            [data_ax_p, row_loc, col_loc, YLim_editboxes_col_loc] = ...
                app.create_new_data_axes_pos(nLabels);
            data_ax_view.make_data_axes(data_ax_p, row_loc, col_loc, data_labels, true(nLabels,1));

            %%%%%% connect all axes to controlers
            data_axes_dict = data_ax_view.data_axes;
            data_axes = data_axes_dict.values;
            data_axes = [data_axes{:}];

            % time_window_controlers
            app.controlers.time_window_controlers.set_controlers(...
                'axes', [data_axes app.views.event_map.local_map_ax]);

            % Y_lim_controlers
            editboxes_in = [num2cell(data_ax_p), row_loc, YLim_editboxes_col_loc];
            editboxes_in = cell_dict(data_labels, num2cell(editboxes_in,2)');
            checkbox_flag = true(size(data_labels));
            app.controlers.Y_lim_controlers.add_ax(data_axes_dict, editboxes_in, checkbox_flag)
            
            % events_controlers
            rows2pos = [row_loc{:}];
            event_ax = app.views.event_map.local_map_ax;
            event_ax_row = event_ax.Layout.Row;
            rows2pos = [rows2pos event_ax_row];
            rows2pos = [min(rows2pos) max(rows2pos)];
            col2pos = unique([col_loc{:}]);
            data_ax_p = unique(data_ax_p);
            app.controlers.events_controlers.set_controlers(...
                "marking_area",{data_ax_p, rows2pos, col2pos})

            % create new events lines to match data
            % use default colors
            % event_types = app.model.events_types;
            % data_ax_view.make_event_lines(event_types);
            % app.views.event_map.update_event_map;

            % make sure event_map start by showing the whole signal
            app.views.event_map.reset_global_map_end();

        end

    end
end