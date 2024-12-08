classdef time_window_controlers_fixed < time_window_controlers
    % same as time_window_controlers, but prevent any change of time
    % windows (only change center)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Setup %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function set_controlers(app, options)
            % same as superclass, but remove the option for win_jumper
            arguments
                app
                
                options.ROI_time
                options.axes

                options.time_jumper
                % options.win_jumper - Removed!
                options.play_speed
                options.play_button
            end
            name_val_opt = namedargs2cell(options);

            set_controlers@time_window_controlers(app, name_val_opt{:})
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Editboxes %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        % Value changed function: time_jumper
        function time_jumper_val_changed(app, evt)
            % same as super, but doesn't let two element  windows

            % stop playing
            app.play_button.Value = 0;
            pause(0.1)

            % collect new jumper value
            user_in = app.time_jumper.Value;

            % parse text 2 numbers
            in_parsed = parse_2num_editbox(user_in);

            % replace infs with maxmimal sig time or zero
            t_end     = app.model.longest_data;
            in_parsed(in_parsed == inf)  = t_end;
            in_parsed(in_parsed == -inf) = 0;

            % validate input
            if numel(in_parsed) > 1 || any(isnan(in_parsed)) || ~issorted(in_parsed,"strictmonotonic")
                % wrong num of elements or bad transfer - flash red, restore
                app.time_jumper.BackgroundColor = 'r';
                pause(0.1)
                app.time_jumper.BackgroundColor = 'w';
                app.time_jumper.Value = evt.PreviousValue;
            else
                % flash green, sort values, change time, make text look nice
                app.time_jumper.BackgroundColor = 'g';
                app.move_window(in_parsed, false)
                pause(0.1)
                app.time_jumper.BackgroundColor = 'w';
                
                % make sure restore view jumps to this point now
                app.set_x_reset_values();

                % move focus to figure
                focus(app.app_base.UIFigure)
            end
        end
    
        % Value changed function: play_button
        function play_button_value_changed(app)
            % zoom out before playing
            if app.play_button.Value ~= 0 % only when starting playing
                zoom(app.app_base.UIFigure, "out")
            end
            play_button_value_changed@time_window_controlers(app)
        end
    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Listners Callbacks %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function update_controls(app)
            % only time_jumper_val_changed when model "in_win_change", and
            % show only new center time
            
            new_time_center = mean(app.model.events_time2disp);
            if ~isempty(app.time_jumper) && isvalid(app.time_jumper)
                app.time_jumper.Value = num2str(new_time_center,'%.3f');
            end
        end
    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Helper Functions %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function move_window(app, new_time, ~)
            % same as superclass, but only let move centers
            
            new_time = mean(new_time);
            
            % deactivate XLims listeners, to prevent multi-calls
            if ~isempty(app.Xlim_changed_listener)
                [app.Xlim_changed_listener.Enabled] = deal(false);
            end

            % move time
            app.model.move_time(new_time)

            % raactivate XLims listeners
            if ~isempty(app.Xlim_changed_listener)
                [app.Xlim_changed_listener.Enabled] = deal(true);
            end
        end

    end
end