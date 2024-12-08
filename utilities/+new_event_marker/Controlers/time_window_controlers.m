classdef time_window_controlers < handle
    %TIME_WINDOW_CONTROLERS - Controllers for managing time-window interactions in a Model-View-Controller (MVC) interface.
    %   This class provides functionality for managing time-related controls such as play button, time jumper,
    %   window jumper, etc., within a MATLAB UI application. 
    %   It serves as part of the controller component in the MVC architectural pattern. 
    %
    %   USAGE:
    %   timeControl = time_movement_controlers(app_base, model, 'Name', Value)
    %
    %   INPUTS:
    %   - app_base: The handle to the parent MATLAB App, holding the UI figure. Requires prop "UIFigure".
    %   - model: An instance of the data2disp_model class representing the model collecting
    %            data & events from sources, and managing what data & events apear in display.
    %
    %   NAME-VALUE PAIR ARGUMENTS:
    %   Those are all specific controlers, that are used to interact with
    %   the time window displayed.
    %   OPTIONAL - skipping any controler will not create it, and will not
    %   add any functionality to existing controlers.
    %
    %   - 'ROI_time':    handle to an ROI object, that its X position will
    %                    define current time window.
    %   - 'axes':        uiaxes handle array, for listening to X-axis limit
    %                    changes & changing current time window accordingly.
    %
    %   The following matlab.ui.controls can be either passed as handles to already existing objects,
    %   or a cell array {[grid layout container parent] [Row] [Col]} specifying
    %   default controler location in base_app.
    %   - 'time_jumper': EditField for manual time adjustment.
    %   - 'win_jumper':  NumericEditField for adjusting the display window size.
    %   - 'play_speed':  NumericEditField for adjusting the playback speed.
    %   - 'play_button': StateButton for starting/stopping playback.
    %
    %   METHODS:
    %   - time_movement_controlers: The constructor method for initializing the time movement controller.
    %   - update_controls: Method for updating controls based on changes in the model.
    %   Callbacks methods:
    %   - time_jumper_val_changed
    %   - win_jumper_val_changed
    %   - play_button_value_changed
    %   - Xlim_changed
    %   - ROI_time_moved
    %
    %   NOTES:
    %   To implement diffrent controlers behaviors, inherite this object
    %   and overwrite callbacks methods.
    %
    %   See also: data2disp_model, add_create_components
    %
    % Author: Lior de Marcas (LdM: https://github.com/Liordemarcas)
    % Date: 2024-02-20


    % adding uicontrols during setup
    % find object to listen to during setup

    % decide controls behavior in advance - overwrite those behavior using
    % inheritance.
    properties
        app_base
        model data2disp_model = data2disp_general_model()
    end

    properties
        controler_axes          matlab.ui.control.UIAxes
        time_jumper             matlab.ui.control.EditField
        win_jumper              matlab.ui.control.NumericEditField
        play_button             matlab.ui.control.StateButton
        play_speed              matlab.ui.control.NumericEditField

        Xlim_changed_listener   event.listener
        ROI_moved_listener      event.listener
        win_changed_listener    event.listener
    end

    properties (Hidden)
        refresh_rate = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Setup %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function app = time_window_controlers(app_base, model)
            % connect controler to model & base

            arguments
                app_base
                model data2disp_model
            end

            % save handles to main objects
            app.app_base = app_base;
            app.model = model;
        end
        
        function set_controlers(app, options)
            arguments
                app
                
                options.ROI_time
                options.axes

                options.time_jumper
                options.win_jumper
                options.play_speed
                options.play_button
            end

            %%%%%%% Connect controlers that do not have defaults

            % add ROI functionality, if needed
            if isfield(options,"ROI_time")
               app.ROI_moved_listener = listener(options.ROI_time,'ROIMoved',@(~,evt) app.ROI_time_moved(evt));
               options = rmfield(options,"ROI_time"); % to skip this object in add_create_components call
            end

            % add reponce to moving x-limits, if needed
            if isfield(options, "axes")
                % create smooth pan responce - note that this must be
                % first, else it reset restore view for some reason
                p = pan(app.app_base.UIFigure);
                p.UseLegacyExplorationModes = "on";
                app.restore_LimitsDimensions_interaction(app.app_base.UIFigure)

                for iAx = 1:numel(options.axes)
                    ax2work = options.axes(iAx);

                    % add listener
                    app.Xlim_changed_listener(iAx) = listener(ax2work,'XLim','PostSet',...
                        @(~,evt) app.Xlim_changed(evt));

                    % set restore view to work smoothly
                    all_buttons = allchild(ax2work.Toolbar);
                    restore_view = all_buttons(strcmpi({all_buttons.Icon}, 'restoreview'));
                    restore_view.ButtonPushedFcn = @(e,d) app.move_on_restore(e,d);
                end

                app.controler_axes = options.axes;
                options = rmfield(options,"axes"); % to skip this object in add_create_components call
            end

            %%%%%%% Connect controlers with defaults

            %%% add Editfields & State Buttons
            % convert options to cell_dict
            options = cell_dict( string(fieldnames(options)), struct2cell(options));

            % create defaults
            components_names = ["time_jumper"; "win_jumper"; "play_speed"; "play_button"];
            def_components = [app.def_control("time_jumper"); app.def_control("win_jumper"); ...
                app.def_control("play_speed"); app.def_control("play_button")];
            
            % create callbacks
            all_callbacks = {@(~,evt) app.time_jumper_val_changed(evt); @(~,~) app.win_jumper_val_changed;...
                ''; @(~,~) app.play_button_value_changed};
            callbacks_names = ["ValueChangedFcn"; "ValueChangedFcn"; ...
                "ValueChangedFcn"; "ValueChangedFcn"];
            
            % create table & push it to "add_create_components"
            components_info_tab = table(def_components, all_callbacks, callbacks_names, ...
                'RowNames',components_names);
            add_create_components(options, app, components_info_tab);

            %%%%%%%% Update controlers on window change in model
            app.win_changed_listener = listener(app.model,"in_win_change", @(~,~) app.update_controls());
        end
    
    end
    
    % helper methods
    methods (Hidden)
        function move_on_restore(~,~,d)
            % make model move when restore view is pressed
            
            % restore view as MATLAB does it
            matlab.graphics.controls.internal.resetHelper(d.Axes,false)

            % make restore trigger axe limit changed
            d.Axes.XLim = d.Axes.XLim;
        end
        
        function restore_LimitsDimensions_interaction(~,fig)
            % recreate axes interaction options "LimitsDimensions"
            % when "UseLegacyExplorationModes" is used
            p = pan(fig);
            z = zoom(fig);
            all_axes = findobj(fig,'Type','axes');

            for iAx = all_axes(:)'
                dim_lim = iAx.InteractionOptions.LimitsDimensions;
                if ~contains(dim_lim, "xy") % assume 2d
                    setAxesPanConstraint(p,iAx,dim_lim)
                    setAxesZoomConstraint(z,iAx,dim_lim)
                end
            end
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Editboxes %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        % Value changed function: play_button
        function play_button_value_changed(app)
            % play the signal.
            % signal will be updated app.frame_rate.Value times a secound,
            % size of update (how much the window is moved) is determined by speed.

            % switch to stop icon, if pressed
            if app.play_button.Value == 1
                app.play_button.Icon = 'slrtStopIcon.png';
            else
                app.play_button.Icon = 'slrtRunIcon.png';
                return
            end

            % collect current window info
            win_size = diff(app.model.events_time2disp, 1, 2);
            current_center = mean(app.model.events_time2disp,2);
            
            % decide in what range refresh rate should be
            min_refresh = 5; % [Hz]
            max_refresh = 20; % [Hz] % based on drawnow limitrate
            if isempty(app.refresh_rate)
                % if no previus rate collected, try for maximal
                app.refresh_rate = max_refresh;
            end

            % collect speed value, transfer it jump size
            req_speed = app.play_speed.Value;
            t_jump = (1/app.refresh_rate)*req_speed;

            % warn if speed is too high
            smallest_window = min(win_size,[],1);
            if smallest_window < abs(t_jump)
                answer = uiconfirm(app.app_base.UIFigure,...
                    sprintf(['Warning:\n'...
                    'Selected speed, and display field may create' ...
                    '"jumps" instead of smooth display.\n'...
                    'Do you want to continue?']),...
                    'Jumps Warning',...
                    'Icon','warning',...
                    'Options',["Yes","No"],...
                    'DefaultOption',"No",...
                    'CancelOption',"No");
                if answer == "No"
                    app.play_button.Icon = 'slrtRunIcon.png';
                    app.play_button.Value = 0;
                    return
                end
            end

            % get inital time mov
            current_t = current_center + t_jump;

            % run play
            while app.play_button.Value ~= 0
                tic;

                % move window
                app.move_window(current_t);

                % wait & let everything refresh
                % pause(1/refresh_interval)
                drawnow

                % break if player reach end/start (note that "move_window"
                % will prevent actually breaking limits)
                next_win = app.model.center2window(current_t, "Events");
                if (next_win(2) >= app.model.longest_data) || (next_win(1) <= 0)
                    break
                end
                
                % check refresh rate - modify it as needed
                elapse_time = toc;
                elapse_time  = elapse_time + 0.01; % pretende app is running slower than it is, so refresh rate will be skewed down
                if elapse_time > (1/app.refresh_rate)
                    app.refresh_rate = max(min_refresh,app.refresh_rate.*0.9);
                elseif elapse_time < (1/app.refresh_rate)
                    app.refresh_rate = min(max_refresh,app.refresh_rate.*1.1);
                end
                
                % wait to make sure refresh rate is limited
                pause((1/app.refresh_rate) - elapse_time);
                
                % change time
                req_speed = app.play_speed.Value;
                t_jump = (1/app.refresh_rate)*req_speed;
                current_t = current_t + t_jump;
            end
            app.play_button.Value = 0;

            % flush memory
            drawnow
            
            % make sure restore view jumps to this point now
            app.set_x_reset_values();

            % restore play icon
            app.play_button.Icon = 'slrtRunIcon.png';
        end

        % Value changed function: time_jumper
        function time_jumper_val_changed(app, evt)
            % validate text in box, if OK change all data time 2 match.
            % note that this function force all axes to have the same
            % limits, but can be change for other cases.

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
            if numel(in_parsed) > 2 || any(isnan(in_parsed)) || ~issorted(in_parsed,"strictmonotonic")
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

        % Value changed function: win_jumper
        function win_jumper_val_changed(app)
            % jump view-window to a sepcific length, while keeping center

            % find new window
            win_size = app.win_jumper.Value;
            t_end = app.model.longest_data; % get recording end
            if isinf(win_size)
                win_size = t_end;
                app.win_jumper.Value = win_size;
            end
            current_center = mean(app.model.events_time2disp,2);
            new_time = current_center + [-win_size win_size]/2;

            % move window
            app.move_window(new_time, false)

            % make sure restore view jumps to this point now
            app.set_x_reset_values();

            % move focus to figure
            focus(app.app_base.UIFigure)
        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Listners Callbacks %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function Xlim_changed(app,evt)
            % collect xlimits, 
            changed_axe = evt.AffectedObject;
            new_lims = changed_axe.XLim;
%             new_center = mean(new_lims);
            % new_lims = evt.NewLimits; 

            % update all
            app.move_window(new_lims, false)
        end
        
        function ROI_time_moved(app, evt)
            % update view when ROI was moved

            % stop playing
            app.play_button.Value = 0;
            pause(0.1)
            
            % collect ROI
            ROI = evt.Source;

            % check if ROI was just moved, or shape changed:
            % if it was only moved, window size is not effected by the ROI
            previous_win = diff(evt.PreviousPosition(:,1));
            current_win  = diff(evt.CurrentPosition(:,1));
            new_time = ROI.Position(:,1)';
            if ~ismembertol(previous_win, current_win)
                % use the time limits taken from ROI pos
                app.move_window(new_time, false)
            else
                % use the center of ROI to choose new pos, keep the window
                % size the same
                app.move_window(new_time, true)
            end
        end
        
        function update_controls(app)
            % update controlers when window was moved in the model.
            % The update is based only on events_time2disp.
            % Note that XLims & ROI_time are not update here, they are
            % update by their view objects.
            
            new_time = app.model.events_time2disp;
            win_size = diff(new_time,1);
            if ~isempty(app.time_jumper) && isvalid(app.time_jumper)
                app.time_jumper.Value = num2str(new_time,'[%.3f %.3f]');
            end
            if ~isempty(app.win_jumper) && isvalid(app.win_jumper)
                app.win_jumper.Value = win_size;
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
                    "space -> press play_button" + newline + ...
                    "rightarrow -> move 1 sec to the right" + newline + ...
                    "leftarrow  -> move 1 sec to the left" + newline;
                return
            end

            % create key responce
            current_center = mean(app.model.events_time2disp);
            switch key_pressed
                case "space"
                    app.play_button.Value = ~app.play_button.Value;
                    app.play_button_value_changed();
                case "rightarrow"
                    app.play_button.Value = 0; % stop playing
                    pause(0.1)
                    app.move_window(current_center + 1)
                case "leftarrow"
                    app.play_button.Value = 0; % stop playing
                    pause(0.1)
                    app.move_window(current_center - 1)
            end
        end
    
    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Helper Functions %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Hidden)

        function move_window(app, new_time, only_centers)
            % move 2 new time window. 
            % Make sure that new_time does not extend over data limits.
            % Note that this format assume 1 time for all data sources &
            % event ax.
            %   INPUTS:
            %       app      - caller obj.
            %       new_time - window center (numeric scalar) or window
            %                  [start stop] time (numeric 2 element row).
            %       only_centers - 
            %                  logical flag. If true, will calculate center
            %                  of new_time & use only it. 
            %                  Default: "true"
            
            % set only_centers default
            if ~exist("only_centers","var")
                only_centers = true;
            end

            % convert time to center if needed
            if only_centers
                new_time = mean(new_time);
            end

            % always pass as start-stop - force all 
            new_time = app.model.center2window(new_time,"Events");
            
            % deal with limits: move center so window will fit
            new_time = app.model.refit_window(new_time);
            
            % deactivate XLims listeners, to prevent multi-calls
            if ~isempty(app.Xlim_changed_listener)
                [app.Xlim_changed_listener.Enabled] = deal(false);
            end

            % move time
            % app.model.move_time_event(new_time, false)
            % app.model.move_time_data(new_time)
            app.model.move_time(new_time)

            % raactivate XLims listeners
            if ~isempty(app.Xlim_changed_listener)
                [app.Xlim_changed_listener.Enabled] = deal(true);
            end
        end
        
        function set_x_reset_values(app)
            % set for each axe it current X limits as its reset values.
            % This help keep "restore_view" behavior consistant with app,
            % instead of sticking to the first time "zoom" was used.
            for iAx = app.controler_axes
                if ~isvalid(iAx)
                    continue
                else
                    iAx.InteractionOptions.RestoredXLimits = iAx.XLim;
                end
            end
        end

        function control_obj = def_control(~, obj_name)
            switch obj_name
                case "time_jumper"
                    control_obj = matlab.ui.control.EditField;
                    control_obj.HorizontalAlignment = 'center';
                    control_obj.Tooltip = {'Time to display. Can be one value (keeping window), range & simple expressions'};
                    control_obj.Placeholder = 'Display Time';
                case "win_jumper"
                    control_obj = matlab.ui.control.NumericEditField;
                    control_obj.Limits = [0 inf];
                    control_obj.LowerLimitInclusive = 'off';
                    control_obj.HorizontalAlignment = 'center';
                    control_obj.Tooltip = {'Time window size (stop-start) to display [s]. Keep view center if possiable, changing view window'};
%                     control_obj.Placeholder = 'Time Window [s]';
                case "play_speed"
                    control_obj = matlab.ui.control.NumericEditField;
                    control_obj.Limits = [-inf inf];
                    control_obj.HorizontalAlignment = 'center';
                    control_obj.Tooltip = {'Play Speed'};
                    %control_obj.Placeholder = 'Play Speed';
                    control_obj.Value = 4;
                    control_obj.ValueDisplayFormat = 'x%.3g';
                case "play_button"
                    control_obj = matlab.ui.control.StateButton;
                    control_obj.Icon = 'slrtRunIcon.png';
                    control_obj.IconAlignment = 'center';
                    control_obj.Text = '';
            end
        end
    end


end