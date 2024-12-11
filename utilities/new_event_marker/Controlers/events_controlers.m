classdef events_controlers < handle

    properties
        app_base
        model data2disp_model = data2disp_general_model()
    end

    properties
        event_chooser_dropdown matlab.ui.control.DropDown
        mark_event_button      matlab.ui.control.StateButton

        event_map event_map 
        % with props:
        %   "event_lines" - will use listener to connect on clicked functionality
        %   "local_map_ax" - will be used as parent to 

        
        time_window_controlers time_window_controlers % controler for time jumps, may contain "play_button" prop

        line_pressed_lst               event.listener
        events_types_changed_lst       event.listener
        marking_area_in_win_change_lst event.listener

        ROI_event    images.roi.Line
        marking_area matlab.ui.control.UIAxes
    end

    properties (Hidden, SetObservable=true)
        marking_area_loc cell = {}
        protected_event_types (:,1) string = string.empty()
    end

    methods
        function app = events_controlers(app_base, model, event_map_view)
            % connect controler to model & base

            arguments
                app_base
                model data2disp_model
                event_map_view
            end

            % save handles to main objects
            app.app_base = app_base;
            app.model = model;
            app.event_map = event_map_view;

        end

        function set_controlers(app, varargin)
      
            %%%%% parse inputs
            p = inputParser;
            p.addParameter("mark_event_button",[])
            p.addParameter("event_chooser_dropdown",[])
            p.addParameter("marking_area",[])

            p.addParameter("listen2line_click",true,@(x) isscalar(x) && (isnumeric(x)||islogical(x)) )
            p.addParameter("protected_event_types",string.empty(),@(x) mustBeText(x))
            p.parse(varargin{:})

            options = p.Results;

            % create a dropdown to choose new event type, if requested
            if ~ismember("mark_event_button", p.UsingDefaults)
                app.create_mark_event_button(options.mark_event_button{:})
            end

            % create a dropdown to choose new event type, if requested
            % (default is always use first existing event type)
            if ~ismember("event_chooser_dropdown" , p.UsingDefaults)
                app.create_event_chooser_dropdown(options.event_chooser_dropdown{:})
                app.events_types_changed_lst = listener(app.model, "events_types_changed", @(~,~) app.update_dropdown_items);
            end

            % create an invisable axe for user to mark events on,
            % if requested (default is only event_map.local_map_ax)
            if ~ismember("marking_area", p.UsingDefaults)
                %app.create_marking_area(options.marking_area{:});
                app.marking_area_loc = options.marking_area;
            end

            % save protected types
            if ~ismember("protected_event_types", p.UsingDefaults)
                app.protected_event_types = options.protected_event_types;
            end

            % add listeners
            if options.listen2line_click
                app.line_pressed_lst = listener(app.event_map,"line_pressed", @(~,evt) app.event_clicked(evt));
            end
        end
    
        function rmv_marking_area(app)
            % delete marking area axe - marking will by default be only in
            % event_map.local_map_ax
            delete(app.marking_area)
        end
    end

    methods (Hidden)
        function create_mark_event_button(app, button_parent, button_row, button_col)
            % create a button for start & stop marking events
            
            button = matlab.ui.control.StateButton();
            button.ValueChangedFcn = @(~,~) app.mark_event_button_val_changed;
            button.Tooltip = {'On press create ROI to mark event.'; 'Press again when finished to collect event info'};
            button.Text = 'Mark Event';
            button.Parent = button_parent;
            button.Layout.Row = button_row;
            button.Layout.Column = button_col;
            app.mark_event_button = button;
            
        end

        function create_event_chooser_dropdown(app, dd_parent, dd_row, dd_col)
            % create a dropdown for choosing what event to mark next

            app.event_chooser_dropdown = uidropdown(dd_parent);
            app.event_chooser_dropdown.Editable = "on";
            app.event_chooser_dropdown.Tooltip = 'Choose which event to use for next marking. Edit to add new event type';
            app.event_chooser_dropdown.Placeholder = 'New Event Type';
            app.event_chooser_dropdown.Layout.Row = dd_row;
            app.event_chooser_dropdown.Layout.Column = dd_col;
            app.event_chooser_dropdown.ValueChangedFcn = @(~,evt) app.event_chooser_val_changed(evt);
        end
    
        function create_marking_area(app, ax_parent, ax_row, ax_col)
           % create an invisible axe to use for marking event on
           
           if ~isempty(app.marking_area) && isvalid(app.marking_area)
               delete(app.marking_area)
           end
           
           % deal with layout parent
           if isa(ax_parent, "matlab.ui.container.GridLayout")
               % is row or col were not given, use the whole grid
               if ~exist('ax_row','var') || isempty(ax_row) || any(isnan(ax_row))
                   ax_row = [1 numel(ax_parent.RowHeight)];
               end
               if ~exist('ax_col','var') || isempty(ax_col) || any(isnan(ax_col))
                   ax_col = [1 numel(ax_parent.ColumnWidth)];
               end
           end
           app.marking_area = uiaxes(ax_parent);
           app.marking_area.Visible = "off";
           app.marking_area.XLim = app.model.events_time2disp;
           % app.marking_area_in_win_change_lst = listener(app.model, "in_win_change", @(~,~) app.model_in_win_change_responce);
           if isa(ax_parent, "matlab.ui.container.GridLayout")
               app.marking_area.Layout.Row = ax_row;
               app.marking_area.Layout.Column = ax_col;
           end
        end
    end


    methods
        function update_dropdown_items(app)
            % prevent event_chooser from being a protected type
            event_types = app.model.events_types;
            event_types(ismember(event_types, app.protected_event_types)) = [];
            if isempty(event_types)
                U = matlab.lang.makeUniqueStrings(["Event"; app.protected_event_types], 1);
                event_types = U(1);
            end
            app.event_chooser_dropdown.Items = event_types;
        end
    end

    methods

        % Value changed function: mark_event_button
        function event_chooser_val_changed(app, evt)
            % prevent event_chooser from being reserved-type
            if ismember(app.event_chooser_dropdown.Value, app.protected_event_types)
                uialert(app.app_base.UIFigure, ...
                    "Types {" + join(app.protected_event_types, ", ") + "} are protected. Use a diffrent name.", ...
                    "Protected Type");
                app.event_chooser_dropdown.BackgroundColor = 'r';
                app.event_chooser_dropdown.Value = evt.PreviousValue;
                pause(0.1)
                app.event_chooser_dropdown.BackgroundColor = [0.96 0.96 0.96];
            end
        end

        % Value changed function: mark_event_button
        function mark_event_button_val_changed(app)
            % mark an event and save its start & end time.

            % stop playing the signal, if possible
            tw_controler = app.time_window_controlers;
            if ~isempty(tw_controler) && isvalid(tw_controler) && ...
                    ~isempty(tw_controler.play_button) && isvalid(tw_controler.play_button)
                tw_controler.play_button.Value = 0;
                pause(0.1)
            end

            % check if to create ROI or collect info, by clicked or relesed
            if app.mark_event_button.Value
                % Pressed: create ROI placment

                % On press, create an ROI for user to mark with, than return
                if ~isempty(app.marking_area_loc)
                    app.create_marking_area(app.marking_area_loc{:});
                    app.ROI_event = drawline(app.marking_area,'DrawingArea','unlimited','Color','k');
                else
                    app.ROI_event = drawline(app.event_map.local_map_ax,'DrawingArea','unlimited','Color','k');
                end
                
            else
                % De-press: Collect ROI result

                % if no line exists, assume user aborted
                if ~isvalid(app.ROI_event)
                    return
                end
                
                % on de-press, collect pos info
                event_t = sort(app.ROI_event.Position(:,1))';
                if ~isempty(app.event_chooser_dropdown) && isvalid(app.event_chooser_dropdown)
                    event_type = string(app.event_chooser_dropdown.Value);
                else
                    event_type = app.model.events_types(1);
                end
                
                % delete ROI and marking area if exists
                delete(app.ROI_event)
                if ~isempty(app.marking_area) && isvalid(app.marking_area)
                    delete(app.marking_area)
                end

                % if both X values are the same - this is an error, or user
                % aborted without marking (as DrawingArea is unlimited). don't save.
                if ismembertol(event_t(1),event_t(2))
                    return
                end

                % push event into model
                app.model.add_event(event_t, event_type)
            end
        end
    
        % Button Down function: event_line
        function event_clicked(app,evt)
            % delete event that was right clicked on

            evt = evt.data.evt;

            % ignore all not - right click events
            if evt.Button ~= 3
                return
            end

            % collect event type - its row name in event_lines table
            clicked_line = evt.Source;
            all_lines = app.event_map.all_lines;
            [match_row,~] = find(ismember(all_lines,clicked_line));
            all_types = app.event_map.lines_types;
            event_type = all_types(match_row);
            
            % if event_type is in a protected type, ignore request
            if ismember(event_type, app.protected_event_types)
                return
            end

            % collect time clicked and find if it is in event
            t_clicked = evt.IntersectionPoint(1);
            relv_event = app.model.point2event(t_clicked, event_type);

            % delete is in event. Confirm before delete.
            if any(relv_event)
                event2del_txt = join(...
                    compose("[%.3f, %.3f]", relv_event(1), relv_event(2)),...
                    ", ");
                answer = uiconfirm(app.app_base.UIFigure,...
                    sprintf(...
                    'Really delete event/s %s: %s',...
                    event_type, event2del_txt...
                    ),...
                    'Delete Event Conformation','Options',["Yes","No"],...
                    'DefaultOption',"No",'CancelOption',"No");
                if answer == "Yes"
                    % delete events
                    app.model.rmv_event(t_clicked, event_type)
                end
            end
        end
    
        % Listener to "in_win_change" event: model
        function model_in_win_change_responce(app)
            % update marking_area limits, if exists
            
            if ~isempty(app.marking_area) && isvalid(app.marking_area)
                app.marking_area.XLim = app.model.events_time2disp;
            end
        end
    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Keypress Responce Callbacks %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function [help_text] = key_responce(app, key_pressed)
            % create keyboard shortcuts

            % create help text
            if nargout > 0
                help_text = ...
                    "shift -> press mark_event_button";
                return
            end

            % create key responce
            switch key_pressed
                case "shift"
                    app.mark_event_button.Value = ~app.mark_event_button.Value;
                    app.mark_event_button_val_changed();
            end
        end
    
    
    end
end