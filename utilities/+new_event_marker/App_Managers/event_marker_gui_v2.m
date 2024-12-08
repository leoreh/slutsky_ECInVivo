classdef event_marker_gui < handle
    % EVENT_MARKER_GUI - A GUI application for marking events in data.
    % This class manages the components of the app, their positions, and assigns controllers to them.
    % It also handles the response to key presses and releases.
    %
    % Properties:
    %   UIFigure - The main figure of the GUI.
    %   containers - Struct holding the different UI containers.
    %   views - Struct holding the different view components.
    %   controlers - Struct holding the different controllers.
    %   model - The data model used by the app.
    %   data_loaded_listener - Listener for the data_loaded event.
    %
    % Methods:
    %   event_marker_gui - Constructor for the event_marker_gui class.
    %   build_app - Initializes the app and sets up the data model and views.
    %   create_app_base - Creates the base app components, such as the figure and containers.
    %   init_app_parts - Initializes the app parts and connects the controllers to the views.
    %   rebuild_app - Rebuilds the app when new data is loaded.
    %   key_release_fcn - Handles key release events.
    %   saveobj - Saves the app state.
    %   loadobj - Loads the app state.
    %   delete_graphics_componants - Deletes all controllers, views, and listeners.
    %   delete - Closes the figure when the object is deleted.
    %
    % Developers:
    %   This is used as a base class for the event_marker_gui class, and it can be used directly.
    %   If you want to customize the app, you can inherit from this class and override the methods you want to change.
    %   Override base_app_build to customize how the app is organized.
    %   Override init_app_parts to customize how the app parts are initialized and connected, independent of data.
    %   Override the rebuild_app method to customize how the app rebuilds itself when data is loaded.
    %   Use composition to add your own controllers and views to the app.
    %   Ensure that your custom controllers and views are properly initialized in the base_app_build method.
    

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
            % EVENT_MARKER_GUI - Constructor for the event_marker_gui class.
            % This method initializes the app by calling the build_app method.
            %
            % Usage:
            %   app = event_marker_gui(varargin)
            %
            % Inputs:
            %   varargin - Optional arguments passed to the build_app method.

            app.build_app(varargin{:});
        end
        
        function build_app(app, varargin)
            % BUILD_APP - Initializes the app and sets up the data model and views.
            % This method is called by the constructor to build the app.
            % It initializes the model, views, and controllers, creates the figure & container, and sets their positions.
            % It also pushes the data into the model, which triggers the data_loaded event and calls the rebuild_app method.
            %
            % Usage:
            %   app.build_app(data_model)
            %
            % Inputs:
            %   data_model - The data model to use for the app, as a first & only argument.
            %                This can be an instance of data2disp_model or any subclass of it.
            %                If not provided, a new instance of data2disp_general_model is created, 
            %                and all inputs are passed to its "add_data" method.
            %
            % Developers:
            %   This method declares what the creation order of the app is.
            %   The order in the base class is: 
            %       1. choose model
            %       2. create the base app (containers)
            %       3. init app parts and connect controlers to views
            %       4. push data into the model, which cause reorganize of the app to fit and display the data.
            %   Note that create_app_base and init_app_parts accepts all app inputs.
            %   by default this is ignored and only effect the model, but you can take advantage of this while overriding.


            %%%%%% Init app parts
            % deal with model
            if nargin == 2 && isa(varargin{1},"data2disp_model")
                app.model = varargin{1};
            else
                % solve weird bug - class do not init new data2disp_general_model upon start
                app.model = data2disp_general_model();
            end
            app.create_app_base(varargin{:}); % create the base app - figure, containers, defualt views & controlers, and set their positions
            app.init_app_parts(varargin{:}); % init model, views, controlers, and make sure that they are connected as needed
            app.data_loaded_listener = listener(app.model, "data_loaded", @(~,~) app.rebuild_app); % refresh app when new data is loaded
            
            % make sure that figure is Visible upon break, even if error
            c = onCleanup(@() set(app.UIFigure, "Visible", "on"));
            
            %%%%%% push data into model -
            % this trigger "data_loaded" event,
            % which cause "rebuild_app" callback
            if nargin == 2 && isa(varargin{1},"data2disp_model")
                % data already in model - just trigger responces to it
                notify(app.model, "data_loaded")
                notify(app.model, "events_types_changed")
                notify(app.model, "in_win_change")
            else
                % actuly push data
                app.model.add_data(varargin{:})
            end
            
            zoom(app.UIFigure, "out") % reset zoom for all axes
        end

        function create_app_base(app, varargin)
            % CREATE_APP_BASE - Creates the base app components, such as the figure and containers.
            %
            % Developers:
            %   This method creates the figure and containers for the app.
            %   It is called by the build_app method to create the base app components.
            %   By default, it creates a figure and two grid layouts, one for the data axes and one for the controls.
            %   Override this method to customize the app's base organization.
            %   This method by default does not require anything else to exist in the app,
            %   but you can take advantage of the varargin to make customisations.

            % create figure
            fig = uifigure('Visible', 'off');
            fig.Color = [1 1 1];
            fig.Position = [100 100 640 480];
            fig.Name = 'Event Marker';
            fig.KeyReleaseFcn = @(~,evt) app.key_release_fcn(evt);
            % UIFigure.KeyPressFcn = @(~,evt) app.key_press_fcn(evt); % leaving here for future use
            app.UIFigure = fig;
            
            % create figure menus
            data_menus = uimenu(app.UIFigure,"Text","&Data");
            app.containers.data_menus = data_menus;

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

        function init_app_parts(app, varargin)
            % INIT_APP_PARTS - Initializes the app parts and connects the controllers to the views.
            %
            % Developers:
            %   This method initializes the model, views, and controllers, and connects the controllers to the views.
            %   It is called by the build_app method to initialize the app parts.
            %   By default, it requires the model to be set up and the containers to be created.
            %   If you want to customize the app parts or their connections, override this method,
            %       and either call the superclass and then customize, or call the parts directly.

            %%%%%% Declare app parts - this is the place to:
            %  1) add new parts
            %  2) pass to each part what other parts it should know about.
            %  3) declare global behavior that is not dependent on data.
            %  4) set all parts to their default positions
            
            % deal with views
            app.views.data_ax_view  = data_ax_view(app, app.model);
            app.views.event_map     = event_map(app, app.model);
            
            % deal with controlers
            app.controlers.time_window_controlers   = time_window_controlers(app, app.model);
            app.controlers.Y_lim_controlers         = Y_lim_controlers(app);
            app.controlers.events_controlers        = events_controlers(app, app.model, app.views.event_map);
            app.controlers.general_controlers       = general_controlers(app, app.model);

            %%%%% Set all parts to their default positions - 
            % this is the place to put app parts in certain locations & containers
            
            % Set View Objects
            app.views.data_ax_view.set_view();
            app.views.event_map.set_view(...
                "local_map_loc", {app.containers.data_pannel, 1, 2}, ...
                "global_map_loc", {app.containers.controls_pannel, 1, [1 5]});
            
            % Set controlers that do not depend on view changes
            app.controlers.general_controlers.set_controlers(...
                "save_button", {app.containers.controls_pannel 2 5}, ...
                "add_close_req", true,...
                "add_data_menu", app.containers.data_menus, ...
                "export_data_menu", app.containers.data_menus, ...
                "save_loc_menu", app.containers.data_menus)
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
                "marking_area",{app.containers.data_pannel, [], 2})
        end
        
        % data_loaded listener for: model
        function rebuild_app(app)
            % REBUILD_APP - Rebuilds the app when new data is loaded.
            % This method is called when the data_loaded event is triggered by the model.
            % It:
            % 1) removes the existing axes, creates new axes to match the data, and connects all axes to the controllers.
            % 2) create & connect Y_lim_controlers to each new axes.
            % 3) connects the time_window_controlers to the new axes.
            % 4) resets the event_map to show the whole signal.
            %
            % Developers:
            %   This method is called by the data_loaded listener to rebuild the app when new data is loaded.
            %   It is called once when the app is first created and then whenever new data is loaded.
            %   By default, it requires the model, controlers and views to be set up - 
            %   only finalizing the app upon data load.
            %   If you want to customize the app's response to new data, override this method.
            %   You can also override the create_new_data_axes_pos method to customize the data_axes positions.
            %   This method is called by the data_loaded listener, which is set up in the build_app method.

            %%%%%% remove existing axes & events
            data_ax_view = app.views.data_ax_view;
            data_ax_view.rmv_data_axes();
            app.views.event_map.rmv_map_lines();
            app.controlers.events_controlers.rmv_marking_area();

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
            % make sure axes resize to match event ax on size change
            % for iAx = 1:numel(data_axes)
            %     addlistener(data_axes(iAx),'SizeChanged',@(~,~) app.rematch_axes(data_labels(iAx)));
            % end
            app.UIFigure.AutoResizeChildren = "off";
            app.UIFigure.SizeChangedFcn = @(varargin) app.resize_fun(data_labels, varargin);
            
            %%%%%% connect all axes to controlers
            data_axes_dict = data_ax_view.data_axes;
            data_axes = data_axes_dict.values;
            data_axes = [data_axes{:}];

            % connect time_window_controlers to new axes to deal with axes manupulations that affect time window 
            % (x-axis limits changes, reset zoom)
            app.controlers.time_window_controlers.set_controlers(...
                'axes', [data_axes app.views.event_map.local_map_ax]);

            % Y_lim_controlers
            editboxes_in = [num2cell(data_ax_p), row_loc, YLim_editboxes_col_loc];
            editboxes_in = cell_dict(data_labels, num2cell(editboxes_in,2)');
            checkbox_flag = true(size(data_labels));
            app.controlers.Y_lim_controlers.add_ax(data_axes_dict, editboxes_in, checkbox_flag)

            % make sure event_map start by showing the whole signal
            app.views.event_map.reset_global_map_end();
        end

    end
 
    %%%%%%%%%%%%% Helper Functions
    methods (Hidden)

        function [data_ax_p, row_loc, col_loc, YLim_editboxes_col_loc] = create_new_data_axes_pos(app,nLabels)
            % CREATE_NEW_DATA_AXES_POS - Determines the positions of new data axes.
            % This method calculates where each new data axis should be placed within the data panel.
            % It ensures that each data axis has its own row and shares a column with the event axis.
            %
            % Usage:
            %   [data_ax_p, row_loc, col_loc, YLim_editboxes_col_loc] = create_new_data_axes_pos(app, nLabels)
            %
            % Inputs:
            %   nLabels - The number of data axes to create.
            %
            % Outputs:
            %   data_ax_p - Array of handles to the parent containers for each data axis.
            %   row_loc - Cell array of row locations for each data axis.
            %   col_loc - Cell array of column locations for each data axis.
            %   YLim_editboxes_col_loc - Cell array of column locations for the YLim edit boxes.
            %
            % Developers:
            %   This method is called by the rebuild_app method to determine the positions of new data axes.
            %   It assumes that the data axes are in the same container (.containers.data_pannel) as the event axis and share a column.
            %   It ensures that each data axis has its own row,
            %   and adjusts the row locations if there is an overlap with the event axis by adding a new row.
            %   Override this method if you need to customize the positions of the data axes.
            %   Ensure that the row and column locations are correctly calculated and returned as cell arrays.
            %   This method expects the event_map view and data_pannel container to be set up before it is called.

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