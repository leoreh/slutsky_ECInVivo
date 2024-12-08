classdef reduct_int_controlers < handle
    % controller for reduct_displayer_IED_disp, to allow for interaction with reduct_displayer

    properties
        app_base reduct_displayer_IED_disp
        model IED_data2disp

        event_map event_map 
        % with props:
        %   "event_lines" - will use listener to connect on clicked functionality
        %   "local_map_ax" - will be used as parent to 
    end

    properties
        return2curr_button matlab.ui.control.Button
        return2prev_button matlab.ui.control.Button
        line_pressed_lst   event.listener
        reduct_disp_closed event.listener
    end


    methods
        function app = reduct_int_controlers(app_base, model, event_map_view)
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

        function set_controlers(app, options)

            arguments
                app

                options.return2curr_button
                options.return2prev_button
            end

            if isfield(options, "return2curr_button")
                app.create_return2curr_button(options.return2curr_button{:})
            end
        
            if isfield(options, "return2curr_button")
                app.create_return2prev_button(options.return2prev_button{:})
            end

            app.line_pressed_lst = listener(app.event_map,"line_pressed", @(~,evt) app.event_clicked(evt));
        end
    end

    methods (Hidden)
        function create_return2curr_button(app, button_parent, button_row, button_col)
            button = matlab.ui.control.Button();
            button.ButtonPushedFcn = @(~,~) app.jump2event("Current");
            button.Tooltip = {'Return to the last choosen detection in the reduct_displayer'};
            button.Text = 'Return 2 Curr';
            button.Parent = button_parent;
            button.Layout.Row = button_row;
            button.Layout.Column = button_col;
            app.return2curr_button = button;
        end

        function create_return2prev_button(app, button_parent, button_row, button_col)
            button = matlab.ui.control.Button();
            button.ButtonPushedFcn = @(~,~) app.jump2event("Previous");
            button.Tooltip = {'Return to the 1 before the last choosen detection in the reduct_displayer'};
            button.Text = 'Return 2 Prev';
            button.Parent = button_parent;
            button.Layout.Row = button_row;
            button.Layout.Column = button_col;
            app.return2prev_button = button;
        end
    end

    methods
        function jump2event(app, jump_type)
            if jump_type == "Current"
                detect2mov = app.app_base.current_pos;
            elseif jump_type == "Previous"
                detect2mov = app.app_base.last_pos;
                if isempty(detect2mov)
                    uialert(app.app_base.UIFigure, "No previous exists!","No previous")
                    return
                end
            end

            detect_t = app.model.collect_detection(detect2mov);
            app.model.move_time(detect_t)
        end
        
        function event_clicked(app,evt)
            
            evt = evt.data.evt;

            % ignore all not - right click events
            if evt.Button ~= 3
                return
            end
            
            % if reduct_app is not open, return without doing anything
            if ~isvalid(app.app_base.reduct_app.UIFigure)
                return
            end
            % collect event type - its row name in event_lines table
            clicked_line = evt.Source;
            all_lines = app.event_map.all_lines;
            [match_row,~] = find(ismember(all_lines,clicked_line));
            all_types = app.event_map.lines_types;
            event_type = all_types(match_row);
            
            % if event_type is in a protected type, ignore request
            if ~ismember(event_type, app.model.temp_event_types)
                return
            end

            % collect time clicked and find if it is in event
            t_clicked = evt.IntersectionPoint(1);
            relv_event = app.model.point2event(t_clicked, event_type);

            % find the closest detection that exist in the model
            relv_event_cent = mean(relv_event);
            detect_samp = app.model.find_closest_detection(relv_event_cent);

            % show on main scatter of reduct_app
            app.app_base.reduct_app.wv_on_other_disp_MenuSelected(detect_samp, [])


        end
    end
end