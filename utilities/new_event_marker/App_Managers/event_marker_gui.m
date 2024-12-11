classdef event_marker_gui < handle
    % hold all of app parts, mange thier positions, assign controlers to
    % the componants of the app.
    properties (Transient)
        UIFigure matlab.ui.Figure

        containers (1,1) struct
        % pannel_holder
        % data_pannel
        % controls_pannel

        views (1,1) struct
        % data_ax_view
        % event_map
        
        controlers (1,1) struct
        % time_window_controlers
        % Ylims_controlers
        % evnt_marking_controlers
        % general_app_controlers

        model data2disp_model = data2disp_general_model()

        data_loaded_listener event.listener

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% App Building %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function app = event_marker_gui(varargin)
            
            app.base_app_build(varargin{:})
            
            % make sure that figure is Visible upon break, even if error
            c = onCleanup(@() set(app.UIFigure, "Visible", "on"));

            % push data into model - 
            % this trigger "data_loaded" event,
            % which cause "rebuild_app" callback
            if nargin == 1 && isa(varargin{1},"data2disp_model")
                % data already in model - just trigger responces to it
                notify(app.model, "data_loaded")
                notify(app.model, "events_types_changed")
                notify(app.model, "in_win_change")
            else
                % actuly push data
                app.model.add_data(varargin{:})
            end
            
            zoom(app.UIFigure, "out")
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
            app.views.event_map     = event_map(app, app.model);

            app.controlers.time_window_controlers   = time_window_controlers(app, app.model);
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
            % app.UIFigure.KeyPressFcn = @(~,evt) app.key_press_fcn(evt);
            app.create_containers()
            
            %%%%%% Set View Objects
            app.views.data_ax_view.set_view();
            app.views.event_map.set_view(...
                "local_map_loc", {app.containers.data_pannel, 1, 2}, ...
                "global_map_loc", {app.containers.controls_pannel, 1, [1 5]});
            
            %%%%%% Set controlers that do not depend on view changes
            data_menus = uimenu(app.UIFigure,"Text","&Data");
            app.controlers.general_controlers.set_controlers(...
                "save_button", {app.containers.controls_pannel 2 5}, ...
                "add_close_req", true,...
                "add_data_menu", data_menus, ...
                "export_data_menu", data_menus, ...
                "save_loc_menu", data_menus)
            app.controlers.time_window_controlers.set_controlers(...
                "time_jumper", {app.containers.controls_pannel 2 1}, ...
                "win_jumper" , {app.containers.controls_pannel 2 2}, ...
                "play_speed",  {app.containers.controls_pannel 2 3}, ...
                "play_button", {app.containers.controls_pannel 2 4},  ...
                "ROI_time", app.views.event_map.ROI_time);
            app.controlers.Y_lim_controlers.set_controlers(...
                "zoom_together_tree", {app.containers.data_pannel 1 1})
            app.controlers.events_controlers.set_controlers(...
                "mark_event_button",{app.containers.controls_pannel 2 6}, ...
                "event_chooser_dropdown",{app.containers.controls_pannel 1 6},...
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
            controls_pannel.ColumnWidth = {'0.3x', '0.3x', '0.3x', '0.3x', '0.3x', ... % Area for time2knob
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
            data_ax_view.make_data_axes(data_ax_p, row_loc, col_loc, data_labels, false(nLabels,1));

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

            % make sure event_map start by showing the whole signal
            app.views.event_map.reset_global_map_end();

        end

    end
 
    %%%%%%%%%%%%% Helper Functions
    methods (Hidden)

        function [data_ax_p, row_loc, col_loc, YLim_editboxes_col_loc] = create_new_data_axes_pos(app,nLabels)
            % create where each data axe should go. Assume that data_axes
            % are in the same container as event_ax, and they share a
            % column. Make sure each axe have their own row.

            % create array of handles to parent containers
            data_ax_p = repelem(app.containers.data_pannel, nLabels, 1);
            
            % choose where to place each axe
            event_ax = app.views.event_map.local_map_ax;
            event_ax_row = event_ax.Layout.Row;
            row_loc = (1:nLabels)';
            overlap_axe = row_loc == event_ax_row;
            if any(overlap_axe)
                % if any data_ax overlaps with event_ax, add a location at
                % the end of the rows.
                nOverlaps = sum(overlap_axe);
                row_loc(overlap_axe,1) = (nLabels+1):(nLabels+nOverlaps);
                % row_loc(overlap_axe) = [];
                % row_loc((end+1):(end+nOverlaps),1) = (nLabels+1):(nLabels+nOverlaps);
            end
            % add rows to parent
            app.containers.data_pannel.RowHeight = repelem({'1x'},1,nLabels+1);
            % transform to cell, as items can span over multiple rows in other cases
            row_loc = num2cell(row_loc);
            
            % place all data_axes in the same col as event_ax
            event_ax_col = event_ax.Layout.Column;
            col_loc = repelem(event_ax_col, nLabels, 1);

            % place YLim editboxes 1 col to the left
            YLim_editboxes_col_loc = col_loc-1;

            % transform to cell, as items can span over multiple cols in other cases
            col_loc = num2cell(col_loc);
            YLim_editboxes_col_loc = num2cell(YLim_editboxes_col_loc);
        end
    
        function resize_fun(app, data_labels, varargin)
            for iAx = 1:numel(data_labels)
                app.rematch_axes(data_labels(iAx));
            end
        end

        % SizeChanged listener: for uiaxes
        function rematch_axes(app, iData)
            % make sure axes_grp match event_ax x-axis location when figure
            % is changed
            data_axes = app.views.data_ax_view.data_axes;
            event_ax = app.views.event_map.local_map_ax;
            data_axes{iData}.InnerPosition([1,3]) = event_ax.InnerPosition([1,3]);
            drawnow
            % pause(0.1)
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
    %%%%%%%%%%%%%%%%%% Save & Load %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        function s = saveobj(app)
            % only save data & current location - use it to recreate app later
            s.model = app.model;
        end
    end

    methods (Static)
        function app = loadobj(s)
            app = event_marker_gui(s.model);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Object Deletion %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function delete_graphics_componants(app)
            % delete all controlers, views & listeners

            for iField = fieldnames(app.controlers)'
                delete(app.controlers.(iField{:}))
            end
            for iField = fieldnames(app.views)'
                delete(app.views.(iField{:}))
            end

            delete(app.data_loaded_listener)
        end

        function delete(app)
            % close figure when object is deleted
            delete(app.UIFigure)
        end
            
    end
end